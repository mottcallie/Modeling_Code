% function for simulating AOTU019 + AOTU025 pursuit drive
%
% INPUTS
% noiseLevel - amplitude of randomized noise added to model
% startPos - starting position for target (e.g., 0 or 180)
% AOTU019synapse - modify AOTU019 synapse (excitatory, inhibitory, slow)
% AOTU019strength - modify AOTU019 synapse strength (# range, 1 for normal)
% AOTUoverlap - modify AOTU019 binocular overlap (normal, medium, low, none)
% DNa02trans - determine if a02 input-output is linear or non-linear (linear, non)
% nTime - number of time points (e.g., 60sec)
%
% OUTPUT
% AOTU019 - visual tuning curve
% AOTU025 - visual tuning curve
% timebase - model time
% visobj_history - trajectory of visual object over time for all repeated runs
%
% Original: RW
% Updated:  MC 02/12/24
% Updated:  MC 04/30/24 added strength multiplyer, performs all runs at once
% Updated:  MC 05/06/24 added linear vs nonlinear option for DNa02
%

function [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model( ...
    noiseLevel,startPos,AOTU019synapse,AOTU019strength,AOTU019overlap,DNa02trans,nTime)
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

%% generate basic visual tuning input curves for AOTU neurons
% visual tuning curves of each cell, over 360 deg azimuthal space

% generate AOTU019 tuning
AOTU019R=cos(f19*deg2rad(visobj_position-c19))';
AOTU019L=cos(f19*deg2rad(visobj_position+c19))';
% modify output range
AOTU019R = ELU(AOTU019R);
AOTU019L = ELU(AOTU019L);
% AOTU019R=max(AOTU019R,0);
% AOTU019L=max(AOTU019L,0);
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

% modify AOTU019 binocular overlap (if needed)
% generate shifted visual tunning curves for AOTU019 with varying degrees
% of bilateral overlap around the midline (0)
switch AOTU019overlap
    case "normal" %normal overlap
        oShift = 0;
    case "medium" %less overlap
        oShift = 10;
    case "low" %even less overlap
        oShift = 20;
    case "none" %no overlap
        oShift = 30;
end
AOTU019R = circshift(AOTU019R,oShift);
AOTU019L = circshift(AOTU019L,-oShift);

% modify AOTU019 strength (if needed) by increasing/decreasing tuning curve gain
AOTU019R = AOTU019R .* AOTU019strength;
AOTU019L = AOTU019L .* AOTU019strength;

% store visual tuning curves for output
AOTU019 = [AOTU019L,AOTU019R];
AOTU025 = [AOTU025L,AOTU025R];


%% generate output curve for DNa02

% set range of possible input values for DNa02
DNa02input = linspace(-1, 2, 1000); 

% set range of possible output values for DNa02
switch DNa02trans
    case "linear"
        % output curve based on linear transformation
        DNa02output = linspace(0,1,1000);
    case "non"
        % output curve based on nonlinear transformation
        DNa02output = ELU(DNa02input);
end

%% run simulated pursuit

% initialize
DNa02R_history = zeros(numRuns,nTime); % DNa02R activity over time
DNa02L_history = zeros(numRuns,nTime); % DNa02L activity over time
rot_vel_history = zeros(numRuns,nTime); % fly's rotational velocity over time

% random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
inputNoise = randn(numRuns,nTime); % drawn from a Gaussian
inputNoise = lowpass(inputNoise',fpass,fs)'; % low-pass filter
inputNoise = noiseLevel*zscore(inputNoise); % rescale

% data storage array for object history
startSide = ones(numRuns,1); 
startSide(2:2:end) = -1; %alternate which side the target starts on
visobj_history = zeros(numRuns,nTime); % position of the visual object over time
visobj_history(:,1:2)=repmat((startPos*startSide),1,2); % starting visual object position, in degrees

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
    % update visual object position (so that positive rotational velocity moves the object leftward)
    visobj_history(:,t)=visobj_history(:,t-1) - rot_vel_history(:,t); 
    visobj_history(:,t)=wrapTo180(visobj_history(:,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
end

end