#!/bin/bash

lumFile=Rootfiles/lum_perrun_dimuon.txt
nevFile=Rootfiles/nev_perday_unix_dimuon.txt

type=mid
#processRunFile=Rootfiles/Run14.run.list
#processRunFile=Rootfiles/all.txt
processRunFile=Rootfiles/Run14.production.$type.list
processRuns=`cat $processRunFile`
outfile=Rootfiles/Run14.production.$type.sf.list
cp blank $outfile

totalLum=0
totalEvt=0
counter=0
totalpre=0
for run in $processRuns; do
    #echo $run
    line=`grep $run $lumFile`
    #echo $line
    if [ ${#line} -gt 0 ]; then
	lum=`echo $line | cut -d " " -f 5`
	prescale=`echo $line | cut -d " " -f 6`
	#echo $lum
	starttime=`echo $line | cut -d " " -f 2`
	#echo $starttime
	line2=`grep -B 1 $starttime $nevFile`
	#echo $line2
	if [ $run -eq 15076101 ]; then
	    nevt=`echo $line2 | cut -d " " -f 2 | cut -d "." -f 1`
	else
	    nevtup=`echo $line2 | cut -d " " -f 4 | cut -d "." -f 1`
	    nevtdown=`echo $line2 | cut -d " " -f 2 | cut -d "." -f 1`
	    nevt=`expr $nevtup - $nevtdown`
	    #echo $nevtup $nevtdown 
	fi
	#echo $nevt
	echo $run $nevt $prescale
	if [ $(bc <<< "$prescale > 1") -eq 1 ]; then
	    totalpre=`expr $totalpre + $nevt`
	fi
	echo $totalpre
	totalEvt=$(($totalEvt + $nevt))
	totalLum=`echo $totalLum + $lum | bc`
	Nmb=`echo $lum*6000000 | bc`
	evt=`echo $nevt*1.0 | bc`
	sf=`echo $Nmb/$evt | bc -l`
	echo $sf >> $outfile
	
    fi
    counter=$(($counter + 1))
    if [ $counter -eq 3000 ]; then
	break
    fi
done

totalNmb=`echo $totalLum*6000000 | bc`
totalEvt=`echo $totalEvt*1.0 | bc`
totalsf=`echo $totalNmb/$totalEvt | bc -l`

echo "Process $counter runs"
echo "Nevents = $totalEvt, L = $totalLum"
echo "MB = $totalNmb"
echo "Scale factor = $totalsf"


