% function for simulating connectome AOTU pursuit drive
% with an oscillatory target
%
% INPUTS
% noiTuning - visual receptive field for all relevant neurons (-180 - 180)
% format as AOTU019, AOTU025, and sum total of all other neurons of interest
% noiseLevel - amplitude of randomized noise added to model
% AOTU019synapse - modify AOTU019 synapse (excitatory, inhibitory, slow)
% nTime - number of time points (e.g., 60sec)
%
% OUTPUT
% visobj_history - trajectory of visual object over time for all repeated runs
% rot_vel_history - rotational velocity output for all repeated runs
%
% 06/20/24 - MC adapted from aotu_steering_oscmodel
% 06/25/24 - MC combine all RFs except 19/25
%

function [timebase,visobj_history,rotvel_history] = aotu_steering_oscmodel2(noiTuning,noiseLevel,AOTU019synapse,nTime)
%% initialize
% load simulation variables
model_settings
loopMode = "open"; %closed or open loop

% set simulation time
timebase = linspace(0,nTime/fs,nTime);

%% generate basic visual tuning input curves for AOTU neurons

% visual tuning curves of each cell, over 360 deg azimuthal space
% tuning curves set based on connectome analyses
AOTU019R = noiTuning.AOTU019;
AOTU019L = flip(AOTU019R);

AOTU025R = noiTuning.AOTU025;
AOTU025L = flip(AOTU025R);

noiSumR = noiTuning.sum;
noiSumL = flip(noiSumR);


%% run simulated pursuit

% initialize
DNa02R_history = zeros(numRuns,nTime); % DNa02R activity over time
DNa02L_history = zeros(numRuns,nTime); % DNa02L activity over time
rotvel_history = zeros(numRuns,nTime); % fly's rotational velocity over time

% random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
rng(13) %set seed
inputNoise = randn(numRuns,nTime); % drawn from a Gaussian
inputNoise = lowpass(inputNoise',fpass,fs)'; % low-pass filter
inputNoise = noiseLevel*zscore(inputNoise); % rescale

% generate triangle wave
wave_partial = linspace(0,wave_amp,1*10); %1/4 of a triangle wave
wave_full = [wave_partial, flip(wave_partial(1:end-1)), -wave_partial(2:end), -flip(wave_partial(2:end-1))]; %alternating triangle wave
wave_fullrep = repmat(wave_full,1,ceil(nTime/length(wave_full))); %repeat for sim run time
wave_fullrep = wave_fullrep(1:nTime); %clip to sim run time

% data storage array for object history
visobj_history = repmat(wave_fullrep,numRuns,1);

% for remaining run duration
for t=3:nTime
    % get previous obj position indices for t-1 and t-2
    [~,p1] = ismember(wrapTo180(round(visobj_history(:,t-1))),visobj_position);
    [~,p2] = ismember(wrapTo180(round(visobj_history(:,t-2))),visobj_position);

    % select steering drive based on AOTU019 synapse modifier
    switch AOTU019synapse
        case "inhibitory"
            current_inputR = AOTU025R(p2) - AOTU019L(p1) + noiSumR(p2);
            current_inputL = AOTU025L(p2) - AOTU019R(p1) + noiSumL(p2);
        case "excitatory"
            current_inputR = AOTU025R(p2) + AOTU019R(p1) + noiSumR(p2);
            current_inputL = AOTU025L(p2) + AOTU019L(p1) + noiSumL(p2);
        case "slowinhibitory"
            current_inputR = AOTU025R(p2) - AOTU019L(p2) + noiSumR(p2);
            current_inputL = AOTU025L(p2) - AOTU019R(p2) + noiSumL(p2);
        case "slowexcitatory"
            current_inputR = AOTU025R(p2) + AOTU019R(p2) + noiSumR(p2);
            current_inputL = AOTU025L(p2) + AOTU019L(p2) + noiSumL(p2);
    end

    % update trajectory
    % determine right vs left steering drives
    dr = interp1(DNa02input,1:length(DNa02input),current_inputR,'nearest'); %right
    DNa02R_history(:,t) = DNa02output(dr); %right
    dl = interp1(DNa02input,1:length(DNa02input),current_inputL,'nearest'); %left
    DNa02L_history(:,t) = DNa02output(dl); %left
    % update steering command based on R-L difference
    rotvel_history(:,t) = k*(DNa02R_history(:,t) - DNa02L_history(:,t)) + inputNoise(:,t);
    % if closed-loop, update visual object position (so that positive rotational velocity moves the object leftward)
    if loopMode=="closed"
        visobj_history(:,t)=visobj_history(:,t) + visobj_history(:,t-1) - rotvel_history(:,t);
        visobj_history(:,t)=wrapTo180(visobj_history(:,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
    end
end


end