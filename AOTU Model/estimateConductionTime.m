% estimateConductionTime calculates conduction time of a neuron
% INPUTS:
%   diameter - diameter of the neuron (in cm)
%   length   - length of the neuron (in cm)
%   rm       - membrane resistance (in Ohm·cm)
%   cm       - membrane capacitance (in F/cm²)
%   ri       - intracellular resistivity (in Ohm·cm)
%
% OUTPUT:
%   conduction_time - estimated conduction time (in seconds)
%
% CREATED: 11/02/2024 - MC
%
function conduction_time = estimateConductionTime(diameter, length, rm, cm, ri)
% Step 1: Calculate the length constant (lambda)
lambda = sqrt((diameter * rm) / (4 * ri));

% Step 2: Calculate the time constant (tau)
tau = rm * cm;

% Step 3: Calculate the conduction velocity (v)
velocity = lambda / tau;

% Step 4: Estimate conduction time (t = L / v)
conduction_time = length / velocity;
end
