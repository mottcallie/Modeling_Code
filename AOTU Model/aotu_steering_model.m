% aotu_steering_model
%
% Simulates the steering behavior of a fly based on AOTU019 and AOTU025 neural input without feedforward control.
% The model generates visual object position, rotational velocity, and inputs from various neurons, 
% and computes the corresponding output based on the AOTU019 synapse type (inhibitory or excitatory).
%
% INPUTS:
%   noiTuning        - Struct containing visual tuning curves for different neurons (AOTU019, AOTU025, sum).
%                      It specifies how these neurons respond to object positions in the azimuthal space.
%   noiseLevel       - Noise level added to the rotational velocity signal.
%   startPos         - Initial object position (azimuthal angle) in degrees.
%   AOTU019synapse   - Type of synapse for AOTU019 (e.g., 'inhibitory', 'excitatory').
%   simDuration      - Duration of the simulation in seconds.
%   k                - Gain factor for computing rotational velocity from neural inputs.
%   alpha            - Alpha parameter for the activation function (ELU) of the DNa02 neuron.
%   runSettings      - Struct containing additional model parameters.
%
% OUTPUTS:
%   timebase         - Time vector for the simulation, from 0 to simDuration, based on the simulation rate (fs).
%   visobj_history   - History of the visual object positions over time for each simulation run.
%   input_history    - History of input values to the DNa02 neurons for each simulation run, for right and left sides.
%   rotvel_history   - History of the fly's rotational velocity (in degrees per second) for each simulation run.
%
% CREATED: 10/30/24 - MC
% UPDATED: 11/06/24 - MC added nonlinearity to AOTU019 and AOTU025
%
function [timebase, visobj_history, input_history, rotvel_history] = aotu_steering_model(...
    noiTuning, noiseLevel, startPos, AOTU019synapse, simDuration, k, alpha, runSettings)

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

%% Generate neural parameters
% Generate delays for each neuron in the model
delay_AOTU019 = round(runSettings.AOTU019_delay / 1000 * fs); % Delay for AOTU019 in samples
delay_Others = round(runSettings.Others_delay / 1000 * fs);   % Delay for AOTU025 and noiSum in samples

% Generate input-output for each neuron
% Visual tuning curves for each neuron type over 360 deg azimuthal space
AOTU019R = abs(noiTuning.AOTU019);  % Right side for AOTU019
AOTU019L = flip(AOTU019R);     % Left side (flipped)
AOTU025R = abs(noiTuning.AOTU025);  % Right side for AOTU025
AOTU025L = flip(AOTU025R);     % Left side (flipped)
groupedR = abs(noiTuning.sum);       % Right side for other neurons
groupedL = flip(groupedR);       % Left side (flipped)
if max(AOTU019R)==0
    AOTU019input = linspace(0, 1, 1000); % Input range
else
    AOTU019input = linspace(0, max(AOTU019R)+0.05, 1000); % Input range
end
AOTU019output = ELU(AOTU019input); % Output range
AOTU025input = linspace(0, max(AOTU025R), 1000); % Input range
AOTU025output = ELU(AOTU025input); % Output range
groupinput = linspace(0, max(groupedR), 1000); % Input range
groupoutput = ELU(groupinput); % Output range

ffNorm = max(AOTU025R)+max(groupedR); %Normalize max output of 025+sum

% Motor input-output for DNa02
DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
DNa02output = adjELU(DNa02input, alpha, runSettings.shift); % Output range for downstream neurons

% Set feedforward component (if in use)
scaleFactor = runSettings.ffscale;
if AOTU019synapse=="feedforward"
    feedforward = 1;
    AOTU019synapse = "inhibitory";
elseif AOTU019synapse=="feedforwardctrl"
    feedforward = 2;
    AOTU019synapse = "inhibitory";
else
    feedforward = 0;
end

%% Run the simulation
%Initialize simulation outputs
input_history = zeros(numRuns, nTime, 2); % Inputs to downstream neurons (right/left)
DNa02R_history = zeros(numRuns, nTime);   % DNa02 right side activity
DNa02L_history = zeros(numRuns, nTime);   % DNa02 left side activity
DNa02RLdiff_history = zeros(numRuns, nTime);   % DNa02 right-left difference activity
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

% For each time step
for t = 3:nTime
    % Compute delayed indices for AOTU019, AOTU025, and noiSum
    t_AOTU019 = max(t - delay_AOTU019, 1); % Apply delay for AOTU019
    t_Others = max(t - delay_Others, 1);   % Apply delay for AOTU025 and noiSum

    % Get previous object positions at the delayed times
    [~, p1_AOTU019] = ismember(wrapTo180(round(visobj_history(:, t_AOTU019))), visobj_position);
    [~, p2_Others] = ismember(wrapTo180(round(visobj_history(:, t_Others))), visobj_position);

    % Compute input-output for AOTU025
    r025i =  interp1(AOTU025input, 1:length(AOTU025input), AOTU025R(p2_Others), 'nearest'); %generate 025 input
    l025i =  interp1(AOTU025input, 1:length(AOTU025input), AOTU025L(p2_Others), 'nearest'); %generate 025 input
    r025o = AOTU025output(r025i)'; %generate 025 output
    l025o = AOTU025output(l025i)'; %generate 025 output
    switch feedforward
        case 0 %no feedforward (standard)
            r019i =  interp1(AOTU019input, 1:length(AOTU019input), AOTU019R(p1_AOTU019), 'nearest'); %generate 019 input
            l019i =  interp1(AOTU019input, 1:length(AOTU019input), AOTU019L(p1_AOTU019), 'nearest'); %generate 019 input
        case 1 %feedforward from others
            feedforward_R = (r025o + groupedR(p2_Others))/ffNorm; %generate feedforward
            feedforward_L = (l025o + groupedL(p2_Others))/ffNorm; %generate feedforward
            AOTU019r_wFF = (AOTU019R(p1_AOTU019)*scaleFactor + feedforward_R)./(scaleFactor+1); %scale visual and feedforward
            AOTU019l_wFF = (AOTU019L(p1_AOTU019)*scaleFactor + feedforward_L)./(scaleFactor+1); %scale visual and feedforward
            r019i =  interp1(AOTU019input, 1:length(AOTU019input), AOTU019r_wFF, 'nearest'); %generate 019 input w/feedforward
            l019i =  interp1(AOTU019input, 1:length(AOTU019input), AOTU019l_wFF, 'nearest'); %generate 019 input w/feedforward
        case 2 %feedforward from self
            feedforward_R = AOTU019output(interp1(AOTU019input, 1:length(AOTU019input), AOTU019R(p2_Others), 'nearest')); %generate feedforward
            feedforward_L = AOTU019output(interp1(AOTU019input, 1:length(AOTU019input), AOTU019L(p2_Others), 'nearest')); %generate feedforward
            AOTU019r_wFF = (AOTU019R(p1_AOTU019)*scaleFactor + feedforward_R')./(scaleFactor+1); %scale visual and feedforward
            AOTU019l_wFF = (AOTU019L(p1_AOTU019)*scaleFactor + feedforward_L')./(scaleFactor+1); %scale visual and feedforward
            r019i =  interp1(AOTU019input, 1:length(AOTU019input), AOTU019r_wFF, 'nearest'); %generate 019 input w/feedforward
            l019i =  interp1(AOTU019input, 1:length(AOTU019input), AOTU019l_wFF, 'nearest'); %generate 019 input w/feedforward
    end
    r019o = AOTU019output(r019i)'; %generate 019 output
    l019o = AOTU019output(l019i)'; %generate 019 output

    % Compute steering drive based on AOTU019 synapse type
    switch AOTU019synapse
        case "inhibitory"
            % AOTU019 inhibitory and contralateral
            current_inputR = r025o - l019o + groupedR(p2_Others);
            current_inputL = l025o - r019o + groupedL(p2_Others);
        case "excitatory"
            % AOTU019 excitatory and ipsilateral
            current_inputR = r025o + r019o + groupedR(p2_Others);
            current_inputL = l025o + l019o + groupedL(p2_Others);
    end
    % Store the input history (right and left inputs)
    input_history(:, t, 1) = current_inputR; % Right input
    input_history(:, t, 2) = current_inputL; % Left input

    %% Apply delay for DNa02 output
    % Compute DNa02 activity without delay
    dr = interp1(DNa02input, 1:length(DNa02input), current_inputR, 'nearest');
    DNa02R_history(:, t) = DNa02output(dr); % Right side DNa02 activity at time t
    dl = interp1(DNa02input, 1:length(DNa02input), current_inputL, 'nearest');
    DNa02L_history(:, t) = DNa02output(dl); % Left side DNa02 activity at time t

    % Calculate difference in DNa02 activity
    DNa02RLdiff_history(:,t) = DNa02R_history(:, t) - DNa02L_history(:, t);
    % Compute the rotational velocity at time t based on scaled DNa02 activity
    rotvel_history(:, t) = k * DNa02RLdiff_history(:,t);

    % Update the visual object position with the delayed rotational velocity influence
    visobj_history(:, t) = visobj_history(:, t) + visobj_history(:, t-1) - rotvel_history(:, t);

    % Add noise to the delayed visual object trajectory
    visobj_history(:, t) = visobj_history(:, t) + inputNoise(:, t);

    % Ensure that the visual object position stays within the range [-180, 180]
    visobj_history(:, t) = wrapTo180(visobj_history(:, t));

end

% Generate timebase for the simulation based on actual sampling rate
timebase = linspace(0, simDuration, nTime); % Simulate from 0 to total simulation duration

% Calculate time step based on sampling rate (fs)
dt = 1 / fs;
rotvel_history = rotvel_history / dt;
