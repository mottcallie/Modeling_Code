% interpolateRFs
% Function to interpolate the RF lookup tables to increase their granularity by 10x.
%
% INPUTS:
%   predicted_RF   - Table containing the visual tuning for .sum, .AOTU019, and .AOTU025.
%   runSettings    - Structure containing run settings including visObjPosition.
%
% OUTPUTS:
%   rf_sum_highres        - Interpolated high-resolution version of RF .sum.
%   rf_AOTU019_highres    - Interpolated high-resolution version of RF .AOTU019.
%   rf_AOTU025_highres    - Interpolated high-resolution version of RF .AOTU025.
%   visObjPosition_highres- High-resolution visual object position (10x longer).
%
function [rf_sum_highres, rf_AOTU019_highres, rf_AOTU025_highres, visObjPosition_highres] = interpolateRFs(predicted_RF, runSettings)
% Extract the original visual object position and RF lookup tables from predicted_RF
visObjPosition = runSettings.visObjPosition;
rf_sum = predicted_RF.sum;
rf_AOTU019 = predicted_RF.AOTU019;
rf_AOTU025 = predicted_RF.AOTU025;

% Define the new high-resolution visual object position for interpolation (10x longer)
visObjPosition_highres = linspace(min(visObjPosition), max(visObjPosition), length(visObjPosition) * 10);

% Interpolate the RF lookup tables to increase granularity (10x longer)
rf_sum_highres = interp1(visObjPosition, rf_sum, visObjPosition_highres, 'pchip');
rf_AOTU019_highres = interp1(visObjPosition, rf_AOTU019, visObjPosition_highres, 'pchip');
rf_AOTU025_highres = interp1(visObjPosition, rf_AOTU025, visObjPosition_highres, 'pchip');
end
