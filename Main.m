clc, clearvars, close all;
L  = 60; % (m)
E  = 25*1e9; % (Pa)
g  = 1529; % (kg/m)
Iz = 6.75; % (m4)
z  = 0.05;
ndof = 6; % number of dof
fs = 500; % (Hz)
amps = 0.01*9.81; % (m/s2)
tf = 60*1; % (s)
addednoise_snr = 5:40; % (dB)
[St,In,Out] = SimpleBeam2D(L,E,g,Iz,z,ndof,fs,amps,tf,addednoise_snr);
