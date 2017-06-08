#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

if [ $# -ne 2 ]
then
  echo "Usage: $0 date enssize"
  exit 1
else
#  rst=$1
  date=$1
  enssize=$2
fi

let Nensp1=enssize+1

cd ${WORKDIR}/ASSIM

count=1
while [ $count -lt ${Nensp1} ]
do  
   count2=`echo 00$count | tail -4c`
   
   cp analysis${count2}.nc ${FORDIR}/shoreface_out.nc
   cp analysis${count2}.nc ${FORDIR}/ARCHIVE/shoreface_analysis${count2}_${date}.nc

let count=count+1
done

exit
