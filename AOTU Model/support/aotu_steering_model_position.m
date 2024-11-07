% aotu_steering_model_position
%
% Computes the steering drive (DNa02R-DNa02L activity) across azimuthal object positions based on neural inputs 
% from AOTU019 and AOTU025. This model evaluates steering behavior without temporal simulation, focusing solely 
% on object position space.
%
% INPUTS:
%   noiTuning        - Struct containing visual tuning curves for AOTU019, AOTU025, and noiSum.
%                      Specifies how each neuron type responds to object positions across azimuthal space.
%   AOTU019synapse   - Type of synapse for AOTU019 ('inhibitory', 'excitatory', or 'feedforward').
%                      Determines how AOTU019 influences the steering drive based on synaptic type.
%   alpha            - Alpha parameter for the activation function (ELU) applied to DNa02 neuron inputs.
%   runSettings      - Struct containing model parameters, including:
%                      `visObjPosition` - Array of azimuthal object positions to evaluate,
%                      `DNa02input` - Input range for DNa02 neurons, and
%                      `DNa02output` - Output range after applying ELU activation.
%
% OUTPUTS:
%   steering_drive   - Array of steering drive (DNa02R-DNa02L activity) across `visObjPosition`.
%                      This array represents the steering input based on the azimuthal position of the visual object.
%
% DESCRIPTION:
%   For each azimuthal position in `visObjPosition`, the function calculates neural inputs based on the selected 
%   synapse type:
%       - 'inhibitory' or 'excitatory' adjusts inputs based on AOTU019’s effects.
%       - 'feedforward' computes inputs with additional scaling, representing a feedforward influence on steering.
%   After calculating DNa02R and DNa02L activities using the ELU activation function, the function computes the 
%   steering drive as the difference in activity between DNa02R and DNa02L. The result, `steering_drive`, describes 
%   how steering responses vary with azimuthal object position.
%
% CREATED: 10/31/24 - MC
%
function [steering_drive] = aotu_steering_model_position(noiTuning, AOTU019synapse, alpha, runSettings)

%% Initialize
% Set default synapse type if not specified
if isempty(AOTU019synapse)
    AOTU019synapse = 'inhibitory';  
end

% Extract object position
visobj_position = runSettings.visObjPosition;

% Initialize arrays for activity and steering drive across visobj_position space
numPositions = length(visobj_position);
DNa02R_activity = zeros(1, numPositions);
DNa02L_activity = zeros(1, numPositions);

%% Generate input-output for each neuron
% Visual tuning curves for each neuron type over 360 deg azimuthal space
AOTU019R = noiTuning.AOTU019;  % Right side for AOTU019
AOTU019L = flip(AOTU019R);     % Left side (flipped)

AOTU025R = noiTuning.AOTU025;  % Right side for AOTU025
AOTU025L = flip(AOTU025R);     % Left side (flipped)

noiSumR = noiTuning.sum;       % Right side for other neurons
noiSumL = flip(noiSumR);       % Left side (flipped)

% Motor input-output for DNa02
DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
DNa02output = adjELU(DNa02input, alpha, runSettings.shift); % Output range for downstream neurons

%% Calculate DNa02 inputs across visobj_position
for i = 1:numPositions
    
    % Compute inputs based on AOTU019 synapse type
    switch AOTU019synapse
        case "inhibitory"
            current_inputR = AOTU025R(i) - AOTU019L(i) + noiSumR(i);
            current_inputL = AOTU025L(i) - AOTU019R(i) + noiSumL(i);
        case "excitatory"
            current_inputR = AOTU025R(i) + AOTU019R(i) + noiSumR(i);
            current_inputL = AOTU025L(i) + AOTU019L(i) + noiSumL(i);
        case "feedforward"
            scaleFactor = 0.5;
            % Compute current drive
            current_driveR = max(0,(AOTU025R(i) - AOTU019L(i) + noiSumR(i))*scaleFactor);
            current_driveL = max(0,(AOTU025L(i) - AOTU019R(i) + noiSumL(i))*scaleFactor);

            % Compute input values based on feedforward scaling
            current_inputR = AOTU025R(i) - (AOTU019L(i) + current_driveL) + noiSumR(i);
            current_inputL = AOTU025L(i) - (AOTU019R(i) + current_driveR) + noiSumL(i);
        case "feedforwardctrl"
            % Compute input values based on only AOTU019 scaling
            scaleFactor = runSettings.ffscale;
            right_feedforward = AOTU019R(i) * scaleFactor;
            left_feedforward = AOTU019L(i) * scaleFactor;
            % AOTU019 receives feedforward input of AOTU019 steering drive
            current_inputR = AOTU025R(i) - (AOTU019L(i) + left_feedforward) + noiSumR(i);
            current_inputL = AOTU025L(i) - (AOTU019R(i) + right_feedforward) + noiSumL(i);
    end
    
    % Calculate DNa02 activity using ELU activation and store in arrays
    DNa02R_activity(i) = interp1(DNa02input, DNa02output, current_inputR, 'nearest');
    DNa02L_activity(i) = interp1(DNa02input, DNa02output, current_inputL, 'nearest');
end

% Calculate steering drive as the difference in DNa02 activity across visobj_position
steering_drive = DNa02R_activity - DNa02L_activity;
%plot(visobj_position, steering_drive)

end
