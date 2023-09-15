#!/bin/bash

currentdir=$1
mniinputfile=$2
currentroi=$3
outputfilename=$4
subjid=$5

out3=".mat"
ext=".txt"

cd $currentdir
fslmaths $mniinputfile -bin bin
fslmaths bin.nii.gz -mul $currentroi mtx
#flirt -in mtx.nii.gz -ref $mniinputfile -out mtx_reg_native -interp nearestneighbour -applyxfm -init trans2.mat
fslmeants -i $mniinputfile -o $outputfilename$ext -m mtx.nii.gz
echo "Timecourse $outputfilename extracted."