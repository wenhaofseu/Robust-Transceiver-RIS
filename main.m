clc;clear;close all;
run Setup.m
run GenerateRISchannel.m
run GenerateRISChannel_Mmax128.m
run GenerateChannelHat.m
run GenerateChannelHat_Mmax128.m
run AO_I.m
run AO_J.m
run MSEverusM.m
run AO_J_SumMSE.m
run AO_J_LowerBound.m
run AO_J_UpperBound.m

cd simulation
run BERverusPt.m
run BERverusPt_NonRobust.m
run BERverusPt_SumMSE.m
cd ..
