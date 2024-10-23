% plotRFsTogether
% Function to plot the RF lookup tables for .sum, .AOTU019, and .AOTU025
% together on a single plot against runSettings.visObjPosition.
% This function also includes flipped versions of each RF.
%
% INPUTS:
%   predicted_RF - Table containing the visual tuning for .sum, .AOTU019, and .AOTU025.
%   runSettings - Structure containing run settings including visObjPosition.
%
% OUTPUTS:
%   This function generates a plot of the RF lookup tables and their flipped 
%   versions on a single graph.
%
function plotRFsTogether(predicted_RF, runSettings)
    % Extract visual object position from runSettings
    visObjPosition = runSettings.visObjPosition;

    % Set colors
    otherColor =   "#77AC30";
    AOTU019Color = "#7E2F8E";
    AOTU025Color = "#0072BD";

    % Plot all RF lookup tables and their flipped versions on a single plot
    hold on;  % Hold to plot all on the same axes

    % Plot .sum and its flipped version
    plot(visObjPosition, predicted_RF.sum, 'Color', otherColor, 'DisplayName', 'Other RF', 'LineWidth', 1);
    plot(visObjPosition, flip(predicted_RF.sum), 'Color', otherColor, 'DisplayName', 'Flipped Other RF', 'LineWidth', 1);

    % Plot .AOTU019 and its flipped version
    plot(visObjPosition, predicted_RF.AOTU019, 'Color', AOTU019Color, 'DisplayName', 'AOTU019 RF', 'LineWidth', 1);
    plot(visObjPosition, flip(predicted_RF.AOTU019), 'Color', AOTU019Color, 'DisplayName', 'Flipped AOTU019 RF', 'LineWidth', 1);

    % Plot .AOTU025 and its flipped version
    plot(visObjPosition, predicted_RF.AOTU025, 'Color', AOTU025Color, 'DisplayName', 'AOTU025 RF', 'LineWidth', 1);
    plot(visObjPosition, flip(predicted_RF.AOTU025), 'Color', AOTU025Color, 'DisplayName', 'Flipped AOTU025 RF', 'LineWidth', 1);

    hold off;

    % Customize the plot
    xlabel('Visual Object Position (deg)');
    ylabel('RF Response');
    title('RF Lookup Tables');
    %legend('show', 'Location', 'northwest');
    grid on;
    xlim([-180 180])
    ylim([0 1.65])
end
