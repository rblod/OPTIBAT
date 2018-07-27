#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh


nb=$1

for i in $(seq -f %02g 1 $nb)
do
\cp ${EXECDIR}/shoreface_out.nc ${OBSDIR}/shoreface_out_${i}.nc 
done
