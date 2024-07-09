% aotu_signalspeed_model
% top level script for modeling a comparison of AOTU019 + AOTU025 signal
% propagation estimates, based on estimates from Gouwens, Wilson 2009
%
% CREATED - MC 05/30/24
%
%% initialize
clear; close all

% range of set cell properties
Rm = 15; %kOhm cm2, membrane resistance
Cm = 1.5; %uF/cm2, membrane capacitance
Ri = 160; %Ohm cm, internal resistance

% estimates for AOTU019/25 length and radius
length_19 = 1585; %um, length
radius_19 = 3; %um, radius

length_25 = 630;
radius_25 = 0.2;

%% estimate signal propagation speeds

% calculate velocity of an unmyelinated axon
Tm = Rm * Cm; % membrane time constant

d_19 = (radius_19*10^-5)*2; %diameter (cm)
d_25 = (radius_25*10^-5)*2;

v_19 = sqrt((radius_19*2) / Tm); %velocity (cm/s)
v_25 = sqrt((radius_25*2) / Tm);

% calculate time for signal to propagate down axon
lengthcm_19 = length_19*10^-5; %cm
lengthcm_25 = length_25*10^-5; %cm

time_19 = (lengthcm_19 / v_19) * 1000; % time, msec
time_25 = (lengthcm_25 / v_25) * 1000; % time, msec