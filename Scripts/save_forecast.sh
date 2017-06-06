#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Usage: $0 date enssize"
  exit 1
else
#  rst=$1
  date=$1
  enssize=$2
fi


WORKDIR='/home/esimon/Rachid/'
FORDIR=${WORKDIR}RUN
CASEDIR='mem'


cd ${FORDIR}

let Nensp1=enssize+1

count=1
while [ $count -lt ${Nensp1} ]
do  
   count2=`echo 00$count | tail -4c`
   
   cd ${FORDIR}/${CASEDIR}${count2}
   
   cp shoreface_out.nc ARCHIVE/shoreface_forecast${count2}_${date}.nc

let count=count+1
done

