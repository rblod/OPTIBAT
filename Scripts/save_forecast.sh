#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

if [ $# -ne 4 ]
then
  echo "Usage: $0 date enssize smooth"
  exit 1
else
#  rst=$1
  date=$1
  enssize=$2
  dateobs=$3
  traj=$4
fi

cd ${FORDIR}

let Nensp1=enssize+1

count=1
while [ $count -lt ${Nensp1} ]
do  
   count2=`echo 00$count | tail -4c`
   
   cd ${FORDIR}/${CASEDIR}${count2}
   [ ! -d  ARCHIVE ] && mkdir ARCHIVE
 
  # cp shoreface_out.nc ARCHIVE/shoreface_forecast${count2}_${date}.nc
  # cp ARCHIVE/shoreface_forecast${count2}_${date}.nc shoreface_out.nc
  
   cp jobout.txt ARCHIVE/jobout_${date}.txt
   if [ $traj == 0 ]
   then 
     # forecast
     cp shoreface_out.nc ARCHIVE/shoreface_forecast${count2}_${date}.nc
   else
     # trajectory for analysis
     cp shoreface_out.nc ARCHIVE/shoreface_trajectory${count2}_${dateobs}.nc
   fi

let count=count+1
done

