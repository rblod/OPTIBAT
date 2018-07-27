#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

\rm -rf ${WORKDIR}
mkdir ${WORKDIR}
cd ${WORKDIR}

for i in `seq ${ID_BEG} ${ENSSIZE}`
do
   mem=`echo 00$i | tail -4c`   
   [ ! -d ${CASEDIR}${mem} ] && mkdir ${CASEDIR}${mem}  
   cd ${WORKDIR}
done

cp ${EXECDIR}/shoreface_out.nc .
python ../PREPRO/genebathy.py ${ENSSIZE} shoreface_out.nc ${ENSSIZE} 0

for i in `seq ${ID_BEG} ${ENSSIZE}`
do
     mem=`echo 00$i | tail -4c`
   
     \mv bathy_${mem}.nc ${CASEDIR}${mem}/shoreface_out.nc
      
     cd ${WORKDIR}
done
      
rm shoreface_out.nc