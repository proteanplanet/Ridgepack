#!/bin/bash
#PBS -q transfer
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -l select=1:ncpus=1
#PBS -A NPSCA07935YF5
#PBS -m e
#PBS -M afrobert@nps.edu

# This PBS script is designed to be run on machines to postprocess RASM or CESM
# output from CICE, organizing output in history files into single variables
# so that they be be processed as timeseries within Ridgepack
#
# Andrew Roberts, Naval Postgraduate School, June 2018 (afrobert@nps.edu)

SYEAR=1980
EYEAR=2009

SOURCEDIR='/p/work1/osinski/archive'

if [ ! -d $WORKDIR/procout ]; then
 mkdir $WORKDIR/procout
fi

cd $WORKDIR/procout

for i in \
R2001dGcaaa01d \
; do

 echo $i

 rasm_cice5_series.bash -c $i -d $SOURCEDIR -s $SYEAR -e $EYEAR >| mon.$i.$EYEAR.txt 2>&1

 rasm_cice5_series.bash -c $i -d $SOURCEDIR -s $SYEAR -e $EYEAR -l >> mon.$i.$EYEAR.txt 2>&1

 rasm_cice5_series.bash -c $i -d $SOURCEDIR -s $SYEAR -e $EYEAR -x >> mon.$i.$EYEAR.txt 2>&1

done

exit 0

