#!/bin/bash
. $(cd $(dirname "$0")/..; pwd)/set_path.sh


echo "Creating ${WORKDIR} and ${CASEDIR}${mem} structure"
echo '--------------------------------------------------'
echo
\rm -rf ${WORKDIR}
mkdir ${WORKDIR}
cd ${WORKDIR}

for i in `seq ${ID_BEG} ${ENSSIZE}`
do
   mem=`echo 00$i | tail -4c`
   
   [ ! -d ${CASEDIR}${mem} ] && mkdir ${CASEDIR}${mem}
done

echo "Cleaning ${ASSIMDIR}"
echo '--------------------------------------------------'
echo
\rm ${ASSIMDIR}/analysis???.nc
\rm ${ASSIMDIR}/forecast???.nc
\rm ${ASSIMDIR}/log_assim_*.txt


echo "Initial perturbation of ${EXECDIR}/shoreface_out.nc"
echo '---------------------------------------------------'
echo
cd ${WORKDIR}
cp ${EXECDIR}/shoreface_out.nc .
python ../PREPRO/genebathy.py ${ENSSIZE} shoreface_out.nc ${ENSSIZE} 0
for i in `seq ${ID_BEG} ${ENSSIZE}`
do
     mem=`echo 00$i | tail -4c`
     \mv bathy_${mem}.nc ${CASEDIR}${mem}/shoreface_out.nc
     cd ${WORKDIR}
done
\rm shoreface_out.nc
