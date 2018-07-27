#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

if [ $# -ne 3 ]
then
  echo "Usage: $0 date_deb enssize date_obs"
  exit 1
else
  date=$1
  enssize=$2
  dateobs=$3
fi

cd ${ROOT_DIR}/ASSIM

[ -f forecast001.nc ] && \rm forecast???.nc

let Nensp1=enssize+1

count=1
while [ $count -lt ${Nensp1} ]
do  
   count2=`echo 00$count | tail -4c`
  
  
   if [ $SMOOTH == 0 ]
   then 
     # link to the current state
     ln -fs ${FORDIR}/${CASEDIR}${count2}/shoreface_out.nc forecast${count2}.nc
   else
     # smoother
     #solution at time_obs
     ln -s ${FORDIR}/${CASEDIR}${count2}/ARCHIVE/shoreface_trajectory${count2}_${dateobs}.nc trajectory${count2}_T001.nc
     
     #solution atbeginning time
     \cp ${FORDIR}/${CASEDIR}${count2}/ARCHIVE/shoreface_forecast${count2}_${date}.nc ${FORDIR}/${CASEDIR}${count2}/shoreface_out.nc
     ln -s ${FORDIR}/${CASEDIR}${count2}/shoreface_out.nc forecast${count2}.nc
   fi  
     
     
   \cp forecast${count2}.nc analysis${count2}.nc
   
let count=count+1
done
