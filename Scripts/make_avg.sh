#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

if [ $# -ne 1 ]
then
  echo "Usage: $0 nsbteps"
  exit 1
else
  dend=$1
fi



dend=$1

for (( i=1;i<=${dend};i++ )); do  ncea -O $(ls ${WORKDIR}/mem*/ARCHIVE/shoreface_trajectory???_${i}.nc) for_avg_${i}.nc; done

for (( i=1;i<=${dend};i++ )); do  ncrcat -O $(ls ${WORKDIR}/mem*/ARCHIVE/shoreface_trajectory???_${i}.nc) for_cat_${i}.nc; cdo timstd for_cat_${i}.nc for_std_${i}.nc; done

for (( i=1;i<=${dend};i++ )); do  ncea -O $(ls ${WORKDIR}/mem*/ARCHIVE/shoreface_analysis???_${i}.nc) ana_avg_${i}.nc; done

for (( i=1;i<=${dend};i++ )); do  ncrcat -O $(ls ${WORKDIR}/mem*/ARCHIVE/shoreface_analysis???_${i}.nc) ana_cat_${i}.nc; cdo timstd ana_cat_${i}.nc ana_std_${i}.nc; done