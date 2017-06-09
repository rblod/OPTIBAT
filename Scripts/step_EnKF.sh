#!/bin/bash

. $(cd $(dirname "$0")/..; pwd)/set_path.sh

if [ $# -ne 3 ]
then
  echo "Usage: $0 ndeb steps step_obs"
  exit 1
else
  ndeb=$1
  steps=$2
  step_obs=$3
fi

cd ${WORKDIR}

##############
# Observations
##############
cd ${SCRIPTDIR}

./build_obs.sh ${step_obs}

############
# Forecast
############

# forecast step
cd ${SCRIPTDIR}
  ./Run_ensemble.sh ${ndeb}

################
# ANALYSIS
###############

# save forecast
cd ${SCRIPTDIR}
let date=${ndeb}+${steps}
if [ $SMOOTH == 0 ]
then 
  traj=1
else
  traj=0
fi   
./save_forecast.sh ${date} ${ENSSIZE} ${step_obs} ${traj}

# link files
 cd ${SCRIPTDIR}
./create_forecast.sh ${date} ${ENSSIZE} ${step_obs}

# EnS
cd ${ASSIMDIR}
[ ! -f EnKF ] && cp ${ROOT_DIR}/EnKF-MPI-Waves/EnKF .
./EnKF enkf.prm > log_assim_${step_obs}.txt


ans=`diff forecast001.nc analysis001.nc`
 if [ -z "${ans}" ] 
then
   echo "There has been no update, we quit!!"
   exit 1;
fi

# prepare and lauch new forecast
# let date=${ndeb}+${steps}
cd ${SCRIPTDIR}
 ./analysis2rst.sh ${date} ${ENSSIZE}

if [ $SMOOTH -ne 0 ]
then 
  ./Run_ensemble.sh ${ndeb}
  ./save_forecast.sh ${date} ${ENSSIZE} ${step_obs} 0
fi   

