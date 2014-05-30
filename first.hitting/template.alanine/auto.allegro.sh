#!/bin/bash

source env.sh

cwd=`pwd`
jobname=${cwd##*/}

rm -f allegro.out allegro.err

echo "#!/bin/bash"				>  submit.sh
echo "#PBS -N $jobname"				>> submit.sh
echo "#PBS -o allegro.out"			>> submit.sh
echo "#PBS -e allegro.err"			>> submit.sh
echo "#PBS -M han.wang@fu-berlin.de"		>> submit.sh
echo "#PBS -l walltime=23:00:00"		>> submit.sh
echo "#PBS -l nodes=1:ppn=$nnodes:AMD"		>> submit.sh
echo "#PBS -l pmem=1000mb"			>> submit.sh

echo "cd $cwd"					>> submit.sh
echo "hostname"					>> submit.sh
echo "date"					>> submit.sh
echo "pwd"					>> submit.sh
echo "./run.sh"					>> submit.sh

qsub submit.sh &> jobid

