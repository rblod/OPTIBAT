ENS_SIZE=100

WORKDIR='/home/esimon/Rachid/RUN/'
CASEDIR='mem'

exec='./wkb.exe'

cd ${WORKDIR}
echo "Run_ensemble.sh: begin"

for i in `seq 1 ${ENS_SIZE}`
do
   mem=`echo 00$i | tail -4c`
   cd ${CASEDIR}${mem}
   
   ${exec} > jobout.txt
   
   cd ${WORKDIR}

done

echo "Run_ensemble.sh: done"
