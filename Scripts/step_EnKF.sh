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
  ./Run_ensemble.sh

################
# ANALYSIS
###############


# save forecast
cd ${SCRIPTDIR}
let date=${ndeb}+${steps}
./save_forecast.sh ${date} ${ENSSIZE}

# link files
 cd ${SCRIPTDIR}
./create_forecast.sh ${ENSSIZE}

# EnS
cd ${ASSIMDIR}
[ ! -f EnKF ] && cp ${ROOT_DIR}/EnKF-MPI-Waves/EnKF .
./EnKF enkf.prm


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
