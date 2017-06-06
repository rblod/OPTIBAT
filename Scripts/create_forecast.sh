#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Usage: $0 enssize"
  exit 1
else
  enssize=$1
fi

WORKDIR='/home/esimon/Rachid/'
FORDIR=${WORKDIR}RUN
CASEDIR='mem'


cd ${WORKDIR}ASSIM

rm forecast???.nc

let Nensp1=enssize+1

count=1
while [ $count -lt ${Nensp1} ]
do  
   count2=`echo 00$count | tail -4c`
  
   ln -s ${FORDIR}/mem${count2}/shoreface_out.nc forecast${count2}.nc
  # cp forecast${count2}.nc analysis${count2}.nc
   
let count=count+1
done
