function [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model_FF_OL(...
    noiTuning, noiseLevel, startPos, AOTU019synapse, simDuration, k, alpha, scaleFactor, runSettings)

%% Initialize
% Set default values if any input is empty
if isempty(noiseLevel)
    noiseLevel = 0;  % Default noise level
end
if isempty(startPos)
    startPos = 0;  % Default starting position
end
if isempty(AOTU019synapse)
    AOTU019synapse = 'inhibitory';  % Default synapse type
end
if isempty(simDuration)
    simDuration = 10;  % Default simulation duration in seconds
end
if isempty(k)
    k = runSettings.k;  % Use k from runSettings if not provided
end
if isempty(alpha)
    alpha = runSettings.alpha;  % Use alpha from runSettings if not provided
end

% Extract model parameters from runSettings
visobj_position = runSettings.visObjPosition; % Object position (azimuthal space)
numRuns = runSettings.numRuns;   % Number of simulation runs
fs = runSettings.fs;             % Simulation update rate (Hz)
fpass = runSettings.fpass;       % Lowpass filter cutoff (Hz)

% Compute the number of time steps based on the simulation duration and fs
nTime = round(simDuration * fs); % Total number of time steps based on duration and update rate

%% Generate delays for each neuron in the model
delay_AOTU019 = round(runSettings.AOTU019_delay / 1000 * fs); % Delay for AOTU019 in samples
delay_Others = round(runSettings.Others_delay / 1000 * fs);   % Delay for AOTU025 and noiSum in samples
delay_DNa02 = round(runSettings.DNa02_delay / 1000 * fs);     % Delay for DNa02 output in samples

%% Generate input-output for each neuron
AOTU019R = noiTuning.AOTU019;  % Right side for AOTU019
AOTU019L = flip(AOTU019R);     % Left side (flipped)
AOTU025R = noiTuning.AOTU025;  % Right side for AOTU025
AOTU025L = flip(AOTU025R);     % Left side (flipped)
noiSumR = noiTuning.sum;       % Right side for other neurons
noiSumL = flip(noiSumR);       % Left side (flipped)

% Motor input-output for DNa02
DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
DNa02output = adjELU(DNa02input, alpha, runSettings.shift); % Output range for downstream neurons

%% Initialize simulation outputs
input_history = zeros(numRuns, nTime, 2); % Inputs to downstream neurons (right/left)
DNa02R_history = zeros(numRuns, nTime);   % DNa02 right side activity
DNa02L_history = zeros(numRuns, nTime);   % DNa02 left side activity
rotvel_history = zeros(numRuns, nTime);   % Fly's rotational velocity over time

% Random noise component for rotational velocity (modeling the influence of noise)
rng(13); % Set seed for reproducibility
inputNoise = randn(numRuns, nTime);       % Random noise drawn from Gaussian distribution
inputNoise = lowpass(inputNoise', fpass, fs)'; % Low-pass filter the noise
inputNoise = noiseLevel * zscore(inputNoise);  % Scale the noise

% Initialize visual object position and history
startSide = ones(numRuns, 1); 
startSide(2:2:end) = -1;  % Alternate which side the object starts on
visobj_history = zeros(numRuns, nTime);  % Visual object position history
visobj_history(:, 1:2) = repmat((startPos * startSide), 1, 2); % Set initial position

%% Generate triangle wave for visual object trajectory
wave_amp = 50; % degrees
wave_partial = linspace(0, wave_amp, 4*10); % 1/4 of a triangle wave
wave_full = [wave_partial, flip(wave_partial(1:end-1)), -wave_partial(2:end), -flip(wave_partial(2:end-1))]; % alternating triangle wave
wave_fullrep = repmat(wave_full, 1, ceil(nTime / length(wave_full))); % repeat for sim run time
wave_fullrep = wave_fullrep(1:nTime); % clip to sim run time

% Assign the repeated triangle wave as the visual object position history
visobj_history = repmat(wave_fullrep, numRuns, 1);

%% Run the simulation for each time step
for t = 3:nTime
    % Compute delayed indices for AOTU019, AOTU025, and noiSum
    t_AOTU019 = max(t - delay_AOTU019, 1); % Apply delay for AOTU019
    t_Others = max(t - delay_Others, 1);   % Apply delay for AOTU025 and noiSum
    t_DNa02 = max(t - delay_DNa02, 1);     % Apply DNa02 delay for motor output

    % Get previous object positions at the delayed times
    [~, p1_AOTU019] = ismember(wrapTo180(round(visobj_history(:, t_AOTU019))), visobj_position);
    [~, p2_Others] = ismember(wrapTo180(round(visobj_history(:, t_Others))), visobj_position);

    % Compute steering drive based on AOTU019 synapse type
    turning_right = rotvel_history(:, t_AOTU019) * scaleFactor;
    turning_right(turning_right < 0) = 0;  % Only positive right turns contribute
    turning_left = -rotvel_history(:, t_AOTU019) * scaleFactor;
    turning_left(turning_left < 0) = 0;    % Only positive left turns contribute

    % Compute input values based on feedforward scaling
    current_inputR = AOTU025R(p2_Others) - (AOTU019L(p1_AOTU019) + turning_left) + noiSumR(p2_Others);
    current_inputL = AOTU025L(p2_Others) - (AOTU019R(p1_AOTU019) + turning_right) + noiSumL(p2_Others);
    
    % Ensure inputs do not exceed limits
    current_inputR(current_inputR < -1) = -1;
    current_inputL(current_inputL < -1) = -1;

    % Store the input history (right and left inputs)
    input_history(:, t, 1) = current_inputR; % Right input
    input_history(:, t, 2) = current_inputL; % Left input

    %% Apply delay for DNa02 output
    dr = interp1(DNa02input, 1:length(DNa02input), current_inputR, 'nearest');
    DNa02R_history(:, t_DNa02) = DNa02output(dr); % Right side DNa02 activity (delayed)
    dl = interp1(DNa02input, 1:length(DNa02input), current_inputL, 'nearest');
    DNa02L_history(:, t_DNa02) = DNa02output(dl); % Left side DNa02 activity (delayed)

    % Compute the rotational velocity as the difference between right and left drives
    rotvel_history(:, t) = k * (DNa02R_history(:, t_DNa02) - DNa02L_history(:, t_DNa02)) + inputNoise(:, t);
end

% Generate timebase for the simulation based on actual sampling rate
timebase = linspace(0, simDuration, nTime); % Simulate from 0 to total simulation duration

end
