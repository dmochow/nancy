% test zeroBadSamples
clear all; close all; clc

X=randn(3,10);
Q=2;
xtent=2;

Xout=zeroBadSamples(X,Q,xtent);

Xout