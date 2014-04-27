#!/bin/bash

submit_command=./auto.hlrn.de

while [ 1 -eq 1 ];
do
    if test ! -f jobid; then
	echo "new submission"
	$submit_command
    else
	checkjob `cat jobid ` &> checkjob.out
	if [ $? -eq 1 ]; then
	    echo "failed at checkjob, no job id, something goes wrong, quit"
	    exit
	fi
	status=`grep State checkjob.out | head -n 1`
	# Idle Running Completed
	if echo $status | grep Completed &> /dev/null; then
	    echo "current job completed"
	    if test -f confout.gro; then
		echo "finished simulation, exit!"
		exit
	    else
		echo "continue simulation, after sleeping 2s"
		sleep 2
		$submit_command
	    fi
	else
	    mydate=`date`
	    echo "checked at $mydate, job is either Idel or Running, should wait, sleep 240s"
	    sleep 240
	fi
    fi
    njobs=`qstat  | grep wangh | wc -l | awk '{print $1}'`
    if [ $njobs -gt 100 ]; then
	echo "$njobs, more than 100, something may be wrong, exit"
	exit
    fi
done
