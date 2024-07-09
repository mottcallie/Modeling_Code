% function for simulating AOTU019 + AOTU025 pursuit drive
%
% INPUTS
% noiseLevel - randomized noise added to target input
% startPos - starting position for target
% AOTU019synapse - modify AOTU019 synapse (excitatory, inhibitory, slow)
% AOTU019strength - modify AOTU019 synapse strength (# range, 1 for normal)
% AOTUoverlap - modify AOTU019 binocular overlap (normal, medium, low, none)
% nTime - number of time points (e.g., 60sec)
%
% OUTPUT
% AOTU019 - visual tuning curve
% AOTU025 - visual tuning curve
% timebase
% visobj_history - trajectory of visual object over time for all repeated runs
%
% Original: RW
% Updated:  MC 02/12/24
% Updated:  MC 04/30/24 added strength multiplyer
%

function [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model( ...
    noiseLevel,startPos,AOTU019synapse,AOTU019strength,AOTU019overlap,nTime)
%% initialize

% visual object positions in azimuthal space ranges from -180 deg to +180 deg
visobj_position = linspace(-180, 180, 361);

% visual tuning variables
f19 = 3; %AOTU019 frequency modifier, increase to reduce tuning width
c19 = 35; %AOTU019 center of tuning curve
f25 = 3; %AOTU025 frequency modifier, increase to reduce tuning width
c25 = 70; %AOTU025 center of tuning curve

% simulation variables
numRuns=10000; % number of times to run the simulation  (ideal is >200)

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
% output curve based on nonlinear transformation

DNa02input = linspace(-1, 2, 1000); % range of possible input values for DNa02
DNa02output = ELU(DNa02input); % range of possible output values for DNa0


%% run simulated pursuit

% initialize simulation
inputNoise = zeros(numRuns,nTime); % random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
startSide = ones(numRuns,1); %alternate which side the target starts on
startSide(2:2:end) = -1;

DNa02R_history = zeros(numRuns,nTime); % DNa02R activity over time
DNa02L_history = zeros(numRuns,nTime); % DNa02L activity over time
rot_vel_history = zeros(numRuns,nTime); % fly's rotational velocity over time
visobj_history = zeros(numRuns,nTime); % position of the visual object over time

% generated repeated runs
for r=1:numRuns
    rng(r+2) % initialize random number generator to produce a new frozen noise snippet
    inputNoise(r,:) = randn(nTime,1); % drawn from a Gaussian
    inputNoise(r,:) = lowpass(inputNoise(r,:),fpass,fs); % low-pass filter
    inputNoise(r,:) = noiseLevel*zscore(inputNoise(r,:)); % rescale
    visobj_history(r,1:2)=startPos*startSide(r); % starting visual object position, in degrees
    
    % for remaining run duration
    for t=3:nTime
        v_tminus1=find(visobj_position==wrapTo180(round(visobj_history(r,t-1)))); % get the index value (v) of the visual object position at time t-1
        v_tminus2=find(visobj_position==wrapTo180(round(visobj_history(r,t-2)))); % get the index value (v) of the visual object position at time t-2
        
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
        d = interp1(DNa02input,1:length(DNa02input),current_inputR,'nearest');
        DNa02R_history(r,t) = DNa02output(d);
        d = interp1(DNa02input,1:length(DNa02input),current_inputL,'nearest');
        DNa02L_history(r,t) = DNa02output(d);
        rot_vel_history(r,t) = k*(DNa02R_history(r,t) - DNa02L_history(r,t)) + inputNoise(r,t); % update steering command;
        visobj_history(r,t)=visobj_history(r,t-1) - rot_vel_history(r,t); % update visual object position (so that positive rotational velocity moves the object leftward)
        visobj_history(r,t)=wrapTo180(visobj_history(r,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
    end

end

end