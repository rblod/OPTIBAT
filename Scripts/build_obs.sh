#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

echo 'Building observations'

if [ $# -ne 1 ]
then
  echo "Usage: $0 date"
  exit 1
else
  date=$1
fi


[ ! -d ${RUNDIR} ] && mkdir ${RUNDIR}
cd ${RUNDIR}

 rm shoreface_obs.nc
 date2=`echo 0$date | tail -3c`
# echo ${OBSDIR}/shoreface_out_${date2}.nc
 ln -s ${OBSDIR}/shoreface_out_${date2}.nc shoreface_obs.nc

[ -f log_obs_${date}.txt ] && \rm log_obs_${date}.txt
touch log_obs_${date}.txt

for obstype in ${OBSTYPES}
do
    echo "   $obstype:"  
    cp infile.data.${obstype} infile.data
    ${ROOT_DIR}/Prep_Routines/prep_obs_waves >> log_obs_${date}.txt
    mv observations.uf observations.uf.${obstype}  
done
rm -f ${ASSIMDIR}/observations.uf
touch ${ASSIMDIR}/observations.uf
for obstype in $OBSTYPES
do
    cat observations.uf.${obstype} >> ${ASSIMDIR}/observations.uf
    \mv  observations-${obstype}.nc ${ASSIMDIR}/.
done
