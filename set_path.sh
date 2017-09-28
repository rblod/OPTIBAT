
#- 
# Path for the repositories, defaults is OPTIBAT
ROOT_DIR=$(cd $(dirname "$0")/..; pwd)

WORKDIR=${ROOT_DIR}/RUN
FORDIR=${ROOT_DIR}/RUN
SCRIPTDIR=${ROOT_DIR}/Scripts
RUNDIR=${ROOT_DIR}/OBS
EXECDIR=${ROOT_DIR}/WKB_MODEL
OBSDIR=${ROOT_DIR}/DATA/
ASSIMDIR=${ROOT_DIR}/ASSIM
CASEDIR='mem'

#
#- information for ssimilation cycle
OBSTYPES="OCG EPB"  # variables to Cuse
exec='wkb.exe'      # executable
ENSSIZE=100         # number of ensemble members
ID_BEG=1
SMOOTH=0



