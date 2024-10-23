% model settings to be held consistent across variations

% object space
visobj_position = linspace(-180, 180, 361); % visual object positions in azimuthal space

% simulation parameters
numRuns=20000; % number of times to run the simulation  (ideal is >20k)
k = 120; % free parameter specifying the maximum change in head direction per time step
fs = 10; % simulation update rate, in Hz
fpass = 2; % lowpass cutoff for random component, in Hz

% open loop parameters
wave_amp = 35; %deg

% downstream process parameters
minIn = -1;
maxIn = 3.4;
DNa02input = linspace(minIn, maxIn, 10000); % set range of possible input values for DNa02
DNa02output = ELU(DNa02input); % output curve based on nonlinear transformation
%figure; plot(DNa02input,DNa02output)