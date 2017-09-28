#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

if [ $# -ne 2 ]
then
  echo "Usage: $0 date_begin date_end"
  exit 1
else
  dbeg=$1
  dend=$2
fi

cd ${SCRIPTDIR}

for k in `seq $dbeg $dend`
do
   echo $k 
   let kp1=k+1
#   kp1=0 
    ./step_EnKF.sh $k 1 $kp1
done
