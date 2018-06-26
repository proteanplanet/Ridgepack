#!/bin/bash

usage(){
cat << EOF

NAME 

$0

SYNOPSIS

rasm_series.bash -f "field(s)" file1 file2 file3 ...
rasm_series.bash -f "field(s)" -d "extra ncks options" file1 file2 file3 ...
rasm_series.bash -f "field(s)" -k file1 file2 file3 ...
rasm_series.bash -h


OPTIONS

-f "field(s)"	Specifiy fields to be extracted (required)

-d "options"	Specifies further options for extracting data from each 
		history file during the extraction process (optional)

-b		Add bounds information from netcdf file to output

-c		Check all input files for fields, rather than first only

-k 		Keep temporary processing files (optional)

-h		Optains this manual page.


DESCRIPTION

This script creates netcdf timeseries from individual RASM history files
from CPL, POP, CICE, WRF, CAM and VIC.  The script specifically handles files with 
naming conventions from RASM/CESM that are time slices.  It extracts individual 
variables from time slice history files, and appends them into a single netcdf file
as a time series.  Certain supporting variables are also extracted, including
coordinate variables.

Fields to be extracted are entered with the -f option as a comma separated
list in quotes e.g. rasm_series.bash -f "uvel,vvel" r29RB1a.cice.h.1990*.nc
Which would create a timeseries netcdf files r29RB1a.cice.h.uvel.vvel.1990.nc
in which only the uvel and vvel variables would be extracted from CICE history
files for all of the 1990 files available.

The NetCDF operator ncks is used to extract the variables from each history file, 
and during this process, it may be required to only extract some part of a field, 
such as a single vertical slab from a VIC or POP variable.  This can be done by
including extra ncks commands in quotes using the -d option, e.g.
rasm_series.bash -f "TEMP" -d "-d z_t,1,1" r29RB1a.pop.h.1991-12-*.nc
which would create a timeseries for the top layer of the ocean for December 1991.
Ths horizontal extent could also be limited:
rasm_series.bash -f "TEMP" -d "-d z_t,1,1 -d nlon,200,500" r29RB1a.pop.h.1991-12-*.nc
An example using WRF is provided here for extracting the lowest level horizonal 
velocity, limited to February and March, 1992::
rasm_series.bash -f "U,V" -d "-d bottom_top,1,1" r30RB1k.wrf.h01.1992-0{2,3}*
This will produce a file called r30RB1k.wrf.h01.U.V.1992.nc
For a complete list of the options available in ncks, use ncks -h.

During construction of the timeseries files, a set of temporary history files
are created containing only the extracted fields requested with supporting 
variables.  If required, these temporary files can be kept at the end of running 
this script, instead of removing them.  If this option is required, use the -k option, 
and the temporary files will be placed in a tar ball and kept. 


VERSION 1.0
Written by Andrew Roberts, March 2012, Naval Postgraduate School
Updated by Andrew Roberts, July 2013, Naval Postgraduate School

EOF
}

FIELDS=
SLICEKEEP=0
NSKSDIMARGS=
EXTRAARGS=
addbounds='no'
checkfields='no'

while getopts "kf:d:h" OPTION
do
     case $OPTION in
         f)
             FIELDS=$OPTARG
             ;;
	 d)
	     NSKSDIMARGS=$OPTARG
	     EXTRAARGS=1
             ;;
         k)
	     SLICEKEEP=1
             ;;
         b)
	     addbounds='yes'
             ;;
         c)
	     checkfields='yes'
             ;;
         h)
             usage
             exit 
             ;;
         ?)
             usage
             exit 1
             ;;
     esac
done


# check that a file has been specified
if [[ -z $FIELDS ]] ; then
 echo 'Missing input fields with -f option';
 echo 'Use -h option for more help';
 exit 1
else
 FIELDLIST=`echo $FIELDS | sed s/\,/\./g`;
fi


# get files to be used after options have been taken
shift $(($OPTIND - 1))
if [ $# -eq 0 ] || [ ! -e $1 ]; then 
     echo 'Missing or incorrect input files';
     echo $1
     echo 'Use -h option for more help';
     exit 1
fi


# find matching times
year1=`basename $1 | cut -f4 -d. | cut -f1 -d-`
month1=`basename $1 | cut -f4 -d. | cut -f2 -d-`
day1=`basename $1 | cut -f4 -d. | cut -f3 -d-`
hour1=`basename $1 | cut -f4 -d. | cut -f4 -d-`


# search for common times, stopping search for year, 
# month, day and hour when differences are found.
for i in $@; do
 if [ ! -z $year1 ]; then 
  year=`basename $i | cut -f4 -d. | cut -f1 -d-`
  if [ $year1 != $year ]; then year1=''; month1=''; day1=''; hour1='' ; fi
 fi
 if [ ! -z $month1 ]; then 
  month=`basename $i | cut -f4 -d. | cut -f2 -d-`
  if [ $month1 != $month ]; then month1=''; day1=''; hour1='' ; fi
 fi
 if [ ! -z $day1 ]; then 
  day=`basename $i | cut -f4 -d. | cut -f3 -d-`
  if [ $day1 != $day ]; then day1=''; hour1=''; fi
 fi
 if [ ! -z $hour1 ]; then 
  hour=`basename $i | cut -f4 -d. | cut -f4 -d-`
  if [ $hour1 != $hour ]; then hour1='' ; fi
 fi
done

# Create timestamp for output file
if [ ! -z $hour1 ]; then 
 timestamp=".$year1-$month1-$day1-$hour1";
elif [ ! -z $day1 ]; then 
 timestamp=".$year1-$month1-$day1";
elif [ ! -z $month1 ]; then 
 timestamp=".$year1-$month1";
elif [ ! -z $year1 ]; then 
 timestamp=".$year1";
else
 timestamp="";
fi

# check netcdf packages are in PATH
packages=(ncks ncrcat ncatted ncdump);
for i in ${packages[@]}; do
 leni=${#i};
 x=`which $i`; lenx=${#x};
 stringstart=`expr $lenx - $leni`;
 y=${x:$stringstart};
 if [ $lenx == 0 ] ; then 
  echo " "
  echo "rasm_series.bash won't run"
  echo "Unable to find $i in your path"; exit 1; 
 fi
done


# get model type
modeltype=`basename $1 | cut -f2 -d. | tr [a-z] [A-Z]`


# get list of float variables in each netcdf file, checking as you go
echo "In `dirname $1`:"
for k in $@; do
 xnameint=`ncdump -h $1  | grep -i -e "int " | grep -e "(" | cut -d"(" -f1 | cut -d" " -f2`
 xnamefloat=`ncdump -h $1  | grep -i -e "float " | grep -e "(" | cut -d"(" -f1 | cut -d" " -f2`
 xnamedouble=`ncdump -h $1  | grep -i -e "double " | grep -e "(" | cut -d"(" -f1 | cut -d" " -f2`
 xname="${xnameint[@]} ${xnamefloat[@]} ${xnamedouble[@]}";
 for i in ${FIELDS//,/ }; do
  included=0;
  for j in $xname; do
   if [ $j == $i ]; then included=1; fi 
  done
  echo -ne "\rChecking $modeltype $i configuration in `basename $k`"
  if [ $included -lt 1 ] ; then 
   echo -ne "\rERROR: $i not found in $k as an int,float or double\n"; exit 1; 
  fi
 done

 # Bypass if every single file needs to be checked instead of just the first
 #if [ $checkfields == 'no' ]; then
 break 
 #fi

done
echo " "

# list descriptor and dimension variables and remove them from variable list
if [ $modeltype == 'CICE' ]; then
 tset="time,tmask,tarea,ANGLE";
elif [ $modeltype == 'POP' ]; then
 tset="time,TAREA,UAREA,ANGLE";
elif [ $modeltype == 'VIC' ]; then
 tset="time,Grid_ID,TMask,longitude,latitude";
elif [ $modeltype == 'CPL' ]; then
# tset="time,doma_area,doma_lon,doma_lat,domo_area,domo_lon,domo_lat,doml_lat,doml_lon,domi_lat,domi_lon";
 tset="time";
elif [ $modeltype == 'WRF' ]; then
 tset="Times,XLAT,XLONG,XLAT_U,XLONG_U,XLAT_V,XLONG_V,ZNU,ZNW";
elif [ $modeltype == 'CAM' ]; then
 tset="time,lat,lon";
else
 echo "Unable to process $modeltype history"; exit 1;
fi

# check descriptors for bounds attributes and add the bounds fields to tset if required
if [ $addbounds == 'yes' ]; then
 ncks -O -c -h -a -q -v $tset $1 rasmseriestemp.nc
 boundsx=`ncdump -h rasmseriestemp.nc |  grep -i "bounds =" | sed -e's/^.*= \"\(.*\)\".*;$/\1,/g'`;
 tset=`echo $boundsx $tset | sed -e 's/ //g'`;
 rm -f rasmseriestemp.nc
fi

# diagnostics of final tset output
echo "Adding $tset"

# create individual netcdf time series for each variable except for time, tlon, tlat
k=0;
for i in $@; do

   k=`expr $k + 1`;
   numb=`printf "%05d\n" $k`;

   name1=`basename $i | cut -f1 -d.`
   name2=`basename $i | cut -f2 -d.`
   name3=`basename $i | cut -f3 -d.`
   name4=`basename $i | cut -f4 -d.`
   newfile=`basename $i`

   if [ $k -gt 1 ]; then
	if [ $name1 != $name1old ]; then echo "File types mismatch: $oldfile, $newfile"; exit 1; fi
	if [ $name2 != $name2old ]; then echo "File types mismatch: $oldfile, $newfile"; exit 1; fi
	if [ $name3 != $name3old ]; then echo "File types mismatch: $oldfile, $newfile"; exit 1; fi
	if [ $name4 == $name4old ]; then echo "Dupicate times: $i, $oldfile"; exit 1; fi
   fi

   name1old=$name1;
   name2old=$name2;
   name3old=$name3;
   name4old=$name4;
   oldfiles=$nfile;

   echo -ne "\rProcessing $name4"

   newfileslice=$name1.tmp.$name2.$name3.$FIELDLIST.$name4.nc

   if [ -e $newfileslice ] ; then rm $newfileslice ; fi

   if [ -z $EXTRAARGS ]; then 
    ncks -O -c -h -a -q -v $FIELDS,$tset $i $newfileslice
   else
    ncks -O -c -h -a -q $NSKSDIMARGS -v $FIELDS,$tset $i $newfileslice || exit 1
   fi

# This section was required for CCSM coupler files, but is no longer needed for CESM
# coupler files, which now include time as the dimension - Andrew Roberts June 2013
#   if [ $modeltype == 'CPL' ]; then # splice in time to CPL files
#    for ij in ${FIELDS//,/ }; do
#        str1="double $ij("
#        str2="double $ij(time, "
#        ncdump $newfileslice | sed -e "s#^.$str1# $str2#" | ncgen -o $newfileslice
#        exit 1
#    done
#   fi

   if [ $SLICEKEEP -gt 0 ] ; then
    ncatted -O -h -a title,global,o,c,"$name1 reduced dataset" $newfileslice
    ncatted -O -h -a source,global,o,c,"RACM: $modeltype" $newfileslice
    ncatted -O -h -a history,global,o,c,"File created `date` by $USER" $newfileslice
    ncatted -O -h -a comment2,global,d,c,'' $newfileslice
    ncatted -O -h -a comment3,global,d,c,'' $newfileslice
   fi

done

newfile=$name1.$name2.$name3.$FIELDLIST$timestamp.nc

if [ -e $newfile ] ; then rm $newfile ; fi

ncrcat -h $name1.tmp.$name2.$name3.$FIELDLIST.*.nc $newfile

echo -ne "\rConcatenated $numb files to $newfile"

ncatted -O -h -a title,global,o,c,"$name1 reduced dataset" $newfile
ncatted -O -h -a source,global,o,c,"RACM: $modeltype" $newfileslice
ncatted -O -h -a history,global,o,c,"File created `date` by $USER" $newfile
ncatted -O -h -a comment2,global,d,c,'' $newfile
ncatted -O -h -a comment3,global,d,c,'' $newfile

if [ $SLICEKEEP -gt 0 ] ; then
 echo "Tarring hyperslabs into $name1.$name2.$name3.$FIELDLIST.tar.gz"
 tar -czf $name1.$name2.$name3.$FIELDS.tar.gz $name1.tmp.$name2.$name3.$FIELDLIST.*.nc
else
 echo ", removing slabs to $name4"
fi

rm $name1.tmp.$name2.$name3.$FIELDLIST.*.nc

exit 0

