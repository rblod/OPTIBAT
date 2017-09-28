#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

\rm -rf ${WORKDIR}
mkdir ${WORKDIR}
cd ${WORKDIR}

cp ../Tools/pert_param .
cp ${ASSIMDIR}/pert_param.in .

for i in `seq ${ID_BEG} ${ENSSIZE}`
do
   mem=`echo 00$i | tail -4c`
   
   [ ! -d ${CASEDIR}${mem} ] && mkdir ${CASEDIR}${mem}

   cp ${EXECDIR}/shoreface_out.nc forecast${mem}.nc
   
   cd ${WORKDIR}
done

./pert_param ${ENSSIZE}

for i in `seq ${ID_BEG} ${ENSSIZE}`
do
     mem=`echo 00$i | tail -4c`
   
     mv forecast${mem}.nc ${CASEDIR}${mem}/shoreface_out.nc
      
     cd ${WORKDIR}
done

\rm    pert_param pert_param.in
      
