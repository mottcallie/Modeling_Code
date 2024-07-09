% function for simulating AOTU019 + AOTU025 pursuit drive
% with an oscillatory target
%
% INPUTS
% noiseLevel - amplitude of randomized noise added to model
% AOTU019synapse - modify AOTU019 synapse (excitatory, inhibitory, slow)
% AOTU019strength - modify AOTU019 synapse strength (# range, 1 for normal)
% AOTU025strength - modify AOTU025 synapse strength (# range, 1 for normal)
% nTime - number of time points (e.g., 60sec)
%
% OUTPUT
% AOTU019 - visual tuning curve
% AOTU025 - visual tuning curve
% timebase
% visobj_history - trajectory of visual object over time for all repeated runs
% rot_vel_history - rotational velocity output for all repeated runs
%
% 05/10/24 - MC adapted from aotu_steering_model
% 05/27/24 - MC addedd open v closed loop options
%

function [AOTU019,AOTU025,timebase,visobj_history,rot_vel_history] = aotu_steering_oscmodel( ...
    noiseLevel,AOTU019synapse,AOTU019strength,AOTU025strength,nTime)
%% initialize

% visual object positions in azimuthal space ranges from -180 deg to +180 deg
visobj_position = linspace(-180, 180, 361);

% visual tuning variables
f19 = 3; %AOTU019 frequency modifier, increase to reduce tuning width
c19 = 35; %AOTU019 center of tuning curve
f25 = 3; %AOTU025 frequency modifier, increase to reduce tuning width
c25 = 70; %AOTU025 center of tuning curve

% simulation variables
numRuns=20000; % number of times to run the simulation  (ideal is >20k)

k = 180; % free parameter specifying the maximum change in head direction per time step
fs = 10; % simulation update rate, in Hz
fpass = 2; % lowpass cutoff for random component, in Hz
timebase = linspace(0,nTime/fs,nTime);

loopMode = "open"; %closed or open loop

%% generate basic visual tuning input curves for AOTU neurons
% visual tuning curves of each cell, over 360 deg azimuthal space

% generate AOTU019 tuning
AOTU019R=cos(f19*deg2rad(visobj_position-c19))';
AOTU019L=cos(f19*deg2rad(visobj_position+c19))';
% modify output range
AOTU019R = ELU(AOTU019R);
AOTU019L = ELU(AOTU019L);
% restrict tuning curve to front field of view
iR = find(AOTU019R == 0);
AOTU019R(1:iR(f19-1)) = 0;
AOTU019R(iR(f19):end) = 0;
iL = find(AOTU019L == 0);
AOTU019L(1:iL(1)) = 0;
AOTU019L(iL(f19-1):end) = 0;

% generate AOTU025 tuning
AOTU025R=cos(f25*deg2rad(visobj_position-c25))';
AOTU025L=cos(f25*deg2rad(visobj_position+c25))';
% modify output range
AOTU025R = ELU(AOTU025R);
AOTU025L = ELU(AOTU025L);
% restrict tuning curve to peripheral field of view
iR = find(AOTU025R == 0);
AOTU025R(1:iR(f25-1)) = 0;
AOTU025R(iR(f25):end) = 0;
iL = find(AOTU025L == 0);
AOTU025L(1:iL(1)) = 0;
AOTU025L(iL(f25-1):end) = 0;

% modify AOTU019 strength (if needed) by increasing/decreasing tuning curve gain
AOTU019R = AOTU019R .* AOTU019strength;
AOTU019L = AOTU019L .* AOTU019strength;
% modify AOTU025 strength (if needed) by increasing/decreasing tuning curve gain
AOTU025R = AOTU025R .* AOTU025strength;
AOTU025L = AOTU025L .* AOTU025strength;

% store visual tuning curves for output
AOTU019 = [AOTU019L,AOTU019R];
AOTU025 = [AOTU025L,AOTU025R];


%% generate output curve for DNa02

% set range of possible input values for DNa02
DNa02input = linspace(-1, 2, 1000); 
% output curve based on nonlinear transformation
DNa02output = ELU(DNa02input);

%% run simulated pursuit

% initialize
DNa02R_history = zeros(numRuns,nTime); % DNa02R activity over time
DNa02L_history = zeros(numRuns,nTime); % DNa02L activity over time
rot_vel_history = zeros(numRuns,nTime); % fly's rotational velocity over time

% random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
inputNoise = randn(numRuns,nTime); % drawn from a Gaussian
inputNoise = lowpass(inputNoise',fpass,fs)'; % low-pass filter
inputNoise = noiseLevel*zscore(inputNoise); % rescale

% generate triangle wave
wave_amp = 35; %deg
wave_partial = linspace(0,wave_amp,1*10); %1/4 of a triangle wave
wave_full = [wave_partial, flip(wave_partial(1:end-1)), -wave_partial(2:end), -flip(wave_partial(2:end-1))]; %alternating triangle wave
wave_fullrep = repmat(wave_full,1,ceil(nTime/length(wave_full))); %repeat for sim run time
wave_fullrep = wave_fullrep(1:nTime); %clip to sim run time

% data storage array for object history
visobj_history = repmat(wave_fullrep,numRuns,1);

% for remaining run duration
for t=3:nTime
    % get previous obj position indices for t-1 and t-2
    [~,v_tminus1] = ismember(wrapTo180(round(visobj_history(:,t-1))),visobj_position);
    [~,v_tminus2] = ismember(wrapTo180(round(visobj_history(:,t-2))),visobj_position);

    % select steering drive based on AOTU019 synapse modifier
    switch AOTU019synapse
        case "inhibitory"
            current_inputR = AOTU025R(v_tminus2) - AOTU019L(v_tminus1);
            current_inputL = AOTU025L(v_tminus2) - AOTU019R(v_tminus1);
        case "excitatory"
            current_inputR = AOTU025R(v_tminus2) + AOTU019R(v_tminus1);
            current_inputL = AOTU025L(v_tminus2) + AOTU019L(v_tminus1);
        case "slow"
            current_inputR = AOTU025R(v_tminus2) - AOTU019L(v_tminus2);
            current_inputL = AOTU025L(v_tminus2) - AOTU019R(v_tminus2);
    end

    % update trajectory
    % determine right vs left steering drives
    dr = interp1(DNa02input,1:length(DNa02input),current_inputR,'nearest'); %right
    DNa02R_history(:,t) = DNa02output(dr); %right
    dl = interp1(DNa02input,1:length(DNa02input),current_inputL,'nearest'); %left
    DNa02L_history(:,t) = DNa02output(dl); %left
    % update steering command based on R-L difference
    rot_vel_history(:,t) = k*(DNa02R_history(:,t) - DNa02L_history(:,t)) + inputNoise(:,t);
    % if closed-loop, update visual object position (so that positive rotational velocity moves the object leftward)
    if loopMode=="closed"
        visobj_history(:,t)=visobj_history(:,t) + visobj_history(:,t-1) - rot_vel_history(:,t);
        visobj_history(:,t)=wrapTo180(visobj_history(:,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
    end
end

end