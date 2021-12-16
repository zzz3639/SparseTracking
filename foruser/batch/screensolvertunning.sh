#!/bin/bash

# submit a sparsetracking solver job tunning all parameters
# Usage: ./run datasetname psfwidth psfdecay boundarysize thint corenum

#matlab=/usr/local/MATLAB/R2018a/bin/matlab
matlab=/home/zzz/MATLAB/R2014b/bin/matlab

foldername=$1
psfwidth=$2
psfdecay=$3
boundarysize=$4
thint=$5
corenum=$6

mkdir ${foldername}all
cd ${foldername}all

foldernameparse=${foldername//// }
read -a foldernamearray <<<$foldernameparse
((arrayend=${#foldernamearray[@]}-1))
jobname=${foldernamearray[arrayend]};

screen -dmS ${jobname}all ${matlab} -nodesktop -singleCompThread -r "load ${foldername}.mat;[ trjans, trjansdetail ] = runSparseTrackingTunningAll(moviethis, movietunning, trjtruetunning, ${psfwidth}, ${psfdecay}, ${boundarysize}, ${thint}, ${corenum}, 'all${jobname}'); save('resfix${jobname}.mat'); exit;"

