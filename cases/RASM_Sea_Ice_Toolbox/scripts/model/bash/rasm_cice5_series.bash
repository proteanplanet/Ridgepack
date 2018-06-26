#!/bin/bash


#####################
# Set up CICE and CPL output variables for averaged monthly output
#####################

helpheader=0
DATALOCATION="$WORKDIR/archive"
TYPE='h'
VARIABLES=(aicen vicen vsnon hi aice hs uvel vvel)
MODTAG='cice'
TAG='ice'
FREQ='monthly'
startyear=1980
endyear=2010
hourly=0
daily=0
x=""
checkdata=1


#####################
# Parse options
#####################

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -h|--help)
    helpheader=1
    ;;
    -d|--directory)
    DATALOCATION="$2";
    shift
    ;;
    -l|--long)
    VARIABLES=(aicen vicen hi aice hs uvel vvel Tsfc uatm vatm sice fswdn flwdn snow rain sst sss uocn vocn fswabs flat fsens flwup evap Tair dsnow meltt melts meltb meltl fsalt fswthru strairx strairy strtltx strtlty strcorx strcory strocnx strocny strintx strinty strength dvidtd dvidtt iage alvl vlvl ardg vrdg apond hpond ipond vsnon)
    ;;
    -x|--coupler)
    TYPE='ha'
    VARIABLES=(a2xavg_Sa_pslv)
    MODTAG='cpl'
    TAG='cpl'
    ;;
    -H|--hourly)
    hourly=1;
    TYPE='h_inst'
    FREQ='hourly'
    ;;
    -D|--daily)
    daily=1;
    FREQ='daily'
    ;;
    -c|--case)
    x="$2"
    shift # past argument
    ;;
    -s|--startyear)
    startyear="$2"
    shift # past argument
    ;;
    -e|--endyear)
    endyear="$2"
    shift # past argument
    ;;
    -b|--bypass)
    checkdata=0
    echo "WARNING: Bypassing data check"
    ;;
    *)
    echo "Unknown option: $1" 
    exit 1
    ;;
esac
shift # past argument or value
done

if [ $daily -eq 1 ] && [ $hourly -eq 1 ]; then
 echo "Can't use -H and -D together"
 exit 1
elif [ $TAG == 'ice' ] && [ $hourly -eq 1 ]; then
 VARIABLES=(divu hi hs aice uvel vvel uatm vatm sig1 sig2 shear a11 a12)
fi

if [ -z $x ] && [ $helpheader -eq 0 ]; then
 echo "Missing case name: Use -c "
 exit 1
fi

#####################
# Print help messages
#####################

if [ $helpheader -eq 1 ]; then
cat <<EOF1
NAME   
    rasm_cice5_series.bash - create timeseries for CICE or CPL output

SYNOPSIS 
    rasm_cice5_series.bash [-h|--help] [-d|--directory] [-l|--long]
                           [-x|--coupler] [-H|--hourly] [-D|--daily]
                           [-c|--case] [-s|--startyear] [-e|--endyear]
                           [-b|--bypass] 

    This script creates a time series from CICE or CPL time slice
    history files.

OPTIONS 

   [-c|--case]      - Specify the case for which the time series will be 
                      created (mandatory). 
   [-d|--directory] - Set the directory of the location of the model case
                      directory (optional). The default is $WORKDIR/archive
   [-l|--long]      - Specifies the full list of CICE5 variables to be turned
                      into a timeseries (optional).
   [-x|--coupler]   - Create a coupler time series instead of CICE series 
                      (optional).
   [-H|--hourly]    - Create time series from hourly data instead of monthly
                      history files (optional).
   [-D|--daily]     - Create time series from daily data instead of monthly
                      history files (optional).
   [-s|--startyear] - Specify the start year of the timeseries (optional). 
                      The default is 1980.
   [-e|--endyear]   - Specify the end year of the timeseries (optional). 
                      The default is 2014. 
   [-b|--bypass]    - Bypass checking that all history files are present
                      before creating timeseries (optional)
   [-h|--help]      - Provides this help page.

EXAMPLES
   
   rasm_cice5_series.bash -c wrfphys_05_rsnw105 -s 1980 -e 2009 -l

   rasm_cice5_series.bash -c wrfphys_05_rsnw105 -s 1996 -e 1996 -H

EOF1

exit 1
fi




#####################
# Check intputs
#####################

if [ -d $DATALOCATION/$x/$TAG/hist ]; then
 echo "Case is: $x"
else
 echo "Unable to find $DATALOCATION/$x/$TAG/hist"
 exit 1
fi


#####################
# create appropriate file list
#####################

echo "Using data located in $DATALOCATION/$x"
filelist=''
if [ $hourly -eq 1 ]; then
 if [ $TAG == 'ice' ] && [ $checkdata -eq 1 ]; then 
  check_rasm_history.bash -d $DATALOCATION -H -i -c $x -s $startyear -e $endyear
  [ $? -eq 1 ] && exit
 elif [ $TAG == 'cpl' ] && [ $checkdata -eq 1 ]; then
  check_rasm_history.bash -d $DATALOCATION -H -x -c $x -s $startyear -e $endyear
  [ $? -eq 1 ] && exit
 fi
 for i in `seq $startyear $endyear`; do
  filelist="$filelist `ls $DATALOCATION/$x/$TAG/hist/$x.$MODTAG.$TYPE.$i-[0-9][0-9]-[0-9][0-9]-[0-9][0-9][0-9][0-9][0-9].nc`"
 done
elif [ $daily -eq 1 ]; then
 if [ $TAG == 'ice' ] && [ $checkdata -eq 1 ]; then 
  check_rasm_history.bash -d $DATALOCATION -D -i -c $x -s $startyear -e $endyear
  [ $? -eq 1 ] && exit
 elif [ $TAG == 'cpl' ] && [ $checkdata -eq 1 ]; then
  check_rasm_history.bash -d $DATALOCATION -D -x -c $x -s $startyear -e $endyear
  [ $? -eq 1 ] && exit
 fi
 for i in `seq $startyear $endyear`; do
  filelist="$filelist `ls $DATALOCATION/$x/$TAG/hist/$x.$MODTAG.$TYPE.$i-[0-9][0-9]-[0-9][0-9].nc`"
 done
else
 if [ $TAG == 'ice' ] && [ $checkdata -eq 1 ]; then 
  check_rasm_history.bash -d $DATALOCATION -i -c $x -s $startyear -e $endyear
  [ $? -eq 1 ] && exit
 elif [ $TAG == 'cpl' ] && [ $checkdata -eq 1 ]; then
  check_rasm_history.bash -d $DATALOCATION -x -c $x -s $startyear -e $endyear
  [ $? -eq 1 ] && exit
 fi
 for i in `seq $startyear $endyear`; do
  filelist="$filelist `ls $DATALOCATION/$x/$TAG/hist/$x.$MODTAG.$TYPE.$i-[0-9][0-9].nc`"
 done
fi


#####################
# move to processing directory
#####################

PROCESSLOCATION=$WORKDIR/processing/$x/$TAG/$FREQ
echo "Processing data in $PROCESSLOCATION"
if [ ! -d "$PROCESSLOCATION" ]; then 
  mkdir -p $PROCESSLOCATION
fi
cd $PROCESSLOCATION
[ $? -eq 1 ] && echo "Can''t move to processing directory" && exit 1
 

#####################
# now create series and pass them to NPS 
#####################

for i in ${VARIABLES[@]}; do
 rasm_series.bash -f $i $filelist
done

exit 

