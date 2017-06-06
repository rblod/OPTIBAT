#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Usage: $0 refdate refhour"
  exit 1
else
#  rst=$1
  year=$1
  month=$2
fi

ARCHDIR=/home/esimon/Selime/Run/Data
WORKDIR=/home/esimon/Selime/Run

cd ${WORKDIR}

rm forecast???.nc

count=1
while [ $count -lt 11 ]
do  
   count2=`echo 00$count | tail -4c`
   count3=`echo 00$count | tail -3c`
   cmonth=`echo 00$month | tail -3c`
   
  
   ln -s ${ARCHDIR}/oops_qg.ens.${count}.${year}T00:00:00Z.P1D${month}.nc forecast${count2}.nc
 
   cp ${ARCHDIR}/oops_qg.ens.${count}.${year}T00:00:00Z.P1D${month}.nc analysis${count2}.nc

let count=count+1
done
