#!/bin/bash

runlist=$1

if [ -z $runlist ]; then
    echo "[w] A run list is needed."
    exit
fi
runs=`cat $runlist`
nruns=`cat $runlist | wc -l`
counter=0
for run in $runs; do
    counter=`expr $counter + 1`
    if [ $counter -lt 10 ]; then
	echo -n "$run, "
    else
	echo "$run,"
	counter=0
    fi
done
echo ""
echo "[i] there are $nruns runs"