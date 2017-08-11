#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

cd ${WORKDIR}
echo "Run_ensemble.sh: begin date"

date=$1
date2=`echo 0$date | tail -3c`

for i in `seq 1 ${ENSSIZE}`
do
   mem=`echo 00$i | tail -4c`
   [ ! -d ${CASEDIR}${mem} ] && mkdir ${CASEDIR}${mem}
   cd ${CASEDIR}${mem}
   ln -sf ${EXECDIR}/wkb.exe .
   ln -sf ${EXECDIR}/namelist .
   ln -sf ${OBSDIR}/bryfile_${date2}.nc shoreface_bry.nc
   cp ${EXECDIR}/shoreface_in.nc .
   
   ./${exec} > jobout.txt
   
   cd ${WORKDIR}

done

echo "Run_ensemble.sh: done"
