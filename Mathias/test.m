clear; clc; close all;
rng('default')
rnds = rand(1,100);
trnd = linspace(0,1,100);
fnc = rnds + trnd;
figure
plot(fnc)
figure
cusum_(fnc)
