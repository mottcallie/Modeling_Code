% plotModelSettings
% Function to plot the RF lookup tables for .sum, .AOTU019, and .AOTU025
% together on a single plot against runSettings.visObjPosition. Also plots
% the ELU transformed values on a separate tile.
%
% INPUTS:
%   predicted_RF - Table containing the visual tuning for .sum, .AOTU019, and .AOTU025.
%   runSettings - Structure containing run settings including visObjPosition.
%
% OUTPUTS:
%   This function generates a plot of the RF lookup tables and their flipped
%   versions on a single graph, and the ELU transformation in a separate plot.

function plotModelSettings(predicted_RF, runSettings)
    % Initialize figure and layout
    figure; set(gcf, 'Position', [100 100 900 400]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');

    % Colors for each RF component
    otherColor = "#77AC30";
    AOTU019Color = "#7E2F8E";
    AOTU025Color = "#0072BD";

    % Plot RF lookup tables and their flipped versions in the first tile
    nexttile;
    hold on;

    % Plot RF .sum and its flipped version
    visObjPosition = runSettings.visObjPosition;
    plot(visObjPosition, predicted_RF.sum, 'Color', otherColor, 'DisplayName', 'Other RF', 'LineWidth', 1);
    plot(visObjPosition, flip(predicted_RF.sum), 'Color', otherColor, 'DisplayName', 'L Other RF', 'LineWidth', 1);

    % Plot RF .AOTU019 and its flipped version
    plot(visObjPosition, predicted_RF.AOTU019, 'Color', AOTU019Color, 'DisplayName', 'AOTU019 RF', 'LineWidth', 1);
    plot(visObjPosition, flip(predicted_RF.AOTU019), 'Color', AOTU019Color, 'DisplayName', 'L AOTU019 RF', 'LineWidth', 1);

    % Plot RF .AOTU025 and its flipped version
    plot(visObjPosition, predicted_RF.AOTU025, 'Color', AOTU025Color, 'DisplayName', 'AOTU025 RF', 'LineWidth', 1);
    plot(visObjPosition, flip(predicted_RF.AOTU025), 'Color', AOTU025Color, 'DisplayName', 'L AOTU025 RF', 'LineWidth', 1);

    % Customize the RF plot
    xlabel('Visual Object Position (deg)');
    ylabel('RF Response');
    title('RF Lookup Tables');
    xlim([-180 180]);
    ylim([0 1.65]);
    legend('show', 'Location', 'northwest');
    grid on;
    hold off;

    % Plot ELU transformation in the second tile
    nexttile;
    hold on;

    % Retrieve and transform input values
    input_values = runSettings.DNa02input;
    shift = runSettings.shift;
    output_values = adjELU(input_values, runSettings.alpha, shift);

    % Plot the ELU-transformed values
    plot(input_values, output_values, 'DisplayName', ['\alpha = ', num2str(runSettings.alpha)], 'LineWidth', 1);

    % Customize the ELU plot
    xlabel('Input values');
    ylabel('Transformed values');
    title('ELU');
    legend('show', 'Location', 'northwest');
    grid on;
    axis tight;
    hold off;
end
