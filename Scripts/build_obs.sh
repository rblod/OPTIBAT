#!/bin/bash

echo 'Building observations'

if [ $# -ne 1 ]
then
  echo "Usage: $0 date"
  exit 1
else
  date=$1
fi


WORKDIR='/home/esimon/Rachid/OBS/'
OBSDIR='DATA/'
OBSTYPES="EPB"

cd ${WORKDIR}

# rm shoreface_obs.nc
# ln -s ${OBSDIR}shoreface_obs_${date}.nc shoreface_obs.nc

for obstype in ${OBSTYPES}
do
    echo "   $obstype:"  
    cp infile.data.${obstype} infile.data
    ./prep_obs_waves
    mv observations.uf observations.uf.${obstype}  
done
rm -f observations.uf
touch observations.uf
for obstype in $OBSTYPES
do
    cat observations.uf.${obstype} >> observations.uf
done
