#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

[ ! -d ${WORKDIR} ] && mkdir ${WORKDIR}
cd ${WORKDIR}

for i in `seq ${ID_BEG} ${ENSSIZE}`
do
   mem=`echo 00$i | tail -4c`
   
   [ ! -d ${CASEDIR}${mem} ] && mkdir ${CASEDIR}${mem}
   cd ${CASEDIR}${mem}
   
   
   # model
#   \cp ${EXECDIR}/${exec} .
   #cp ${EXECDIR}namelist .
#   ln -fs ${EXECDIR}/namelist .
   
   
#   cp ${EXECDIR}input_param.txt .
# change here for rel obs path
#   cp   ${OBSDIR}/shoreface_out_00.nc shoreface_out.nc
   \cp   ../../WKB_MODEL/shoreface_out.nc .
   
#   cp ARCHIVE/shoreface_forecast${mem}_1.nc shoreface_out.nc
   
#   mkdir ARCHIVE
   
   cd ${WORKDIR}
done
