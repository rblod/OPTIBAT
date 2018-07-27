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

#cd ${WORKDIR} # useless ???
#${SCRIPTDIR}/create_config_ensemble.sh

if [ "$1" == 1 ]; then
#  [ -f ${SCRIPTDIR}/step_EnKF.output ] && \rm ${SCRIPTDIR}/step_EnKF.output
  echo '--------------------------' #> ${SCRIPTDIR}/step_EnKF.output
  echo '--- Welcome to OPTIBAT ---' ##>> ${SCRIPTDIR}/step_EnKF.output
  echo '--------------------------' ##>> ${SCRIPTDIR}/step_EnKF.output
  echo '  ' ##>> ${SCRIPTDIR}/step_EnKF.output
fi

echo '--- STARTING ASSIMILATION CYCLE  No ---' $ndeb ##>> ${SCRIPTDIR}/step_EnKF.output
echo '  ' #>> ${SCRIPTDIR}/step_EnKF.output


##############
# Observations
##############

cd ${SCRIPTDIR}
echo '--- Building obs  ---' #>> ${SCRIPTDIR}/step_EnKF.output
echo '  ' #>> ${SCRIPTDIR}/step_EnKF.output
./build_obs.sh ${step_obs} #>> step_EnKF.output

############
# Forecast
############

# forecast step
cd ${SCRIPTDIR}
echo '--- Running forecast  ---' #>> ${SCRIPTDIR}/step_EnKF.output
echo '  ' #>> ${SCRIPTDIR}/step_EnKF.output
./Run_ensemble.sh ${ndeb} #>> step_EnKF.output

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
echo '  ' #>> ${SCRIPTDIR}/step_EnKF.output
echo '--- Saving forecast  ---' #>> ${SCRIPTDIR}/step_EnKF.output
echo '  ' #>> ${SCRIPTDIR}/step_EnKF.output
./save_forecast.sh ${date} ${ENSSIZE} ${step_obs} ${traj} #>> step_EnKF.output

# link files
cd ${SCRIPTDIR}
echo '--- Link forecast  ---' #>> ${SCRIPTDIR}/step_EnKF.output
echo "  " #>> ${SCRIPTDIR}/step_EnKF.output
./create_forecast.sh ${date} ${ENSSIZE} ${step_obs} #>> step_EnKF.output

# EnS
cd ${ASSIMDIR}
echo '--- Starting assimilation forecast  ---' #>> ${SCRIPTDIR}/step_EnKF.output
echo '    ' #>> ${SCRIPTDIR}/step_EnKF.output
\cp ${ROOT_DIR}/EnKF-MPI-Waves/EnKF .
./EnKF enkf.prm > log_assim_${step_obs}.txt

ans=`diff forecast001.nc analysis001.nc`
if [ -z "${ans}" ] 
then
   echo 'There has been no update, we quit!!'
   exit 1;
else
   echo 'OK, We start another iteration'   
fi

# prepare and lauch new forecast
# let date=${ndeb}+${steps}
echo '--- Launching new forecast  ---' #>> ${SCRIPTDIR}/step_EnKF.output
echo '  ' #>> ${SCRIPTDIR}/step_EnKF.output
cd ${SCRIPTDIR}
 ./analysis2rst.sh ${date} ${ENSSIZE} #>> step_EnKF.output

if [ $SMOOTH -ne 0 ]
then 
  ./Run_ensemble.sh ${ndeb} #>> step_EnKF.output
  ./save_forecast.sh ${date} ${ENSSIZE} ${step_obs} 0 #>> step_EnKF.output
fi   

