% generateSteeringPlot
%
% Generates a plot of the steering drive across azimuthal object positions for two different conditions.
% The function saves the plot as both PNG and SVG formats.
%
% INPUTS:
%   steering_drive   - 2-row array, output of aotu_steering_model_position. Each row contains the steering drive
%                      data for a different condition.
%   comparisonLabel  - 1x2 cell array containing labels for each of the two conditions (e.g., {'Condition 1', 'Condition 2'}).
%   runSettings      - Struct containing model parameters, including:
%                      `visObjPosition` - Array of azimuthal object positions to plot on the x-axis.
%   folder           - Struct with fields specifying directories for saving final PNG and vector graphics (SVG) files.
%
% DESCRIPTION:
%   This function creates a new figure, plots the steering drive across object positions for two conditions, and
%   saves the plot as PNG and SVG files in the specified directories.
%
% CREATED: 10/31/24 - MC
%
function generateSteeringPlot(steering_drive, comparisonLabel, runSettings, folder)

    % Extract object positions from runSettings
    visObjPosition = runSettings.visObjPosition;
    
    % Generate new figure
    figure;
    hold on;
    
    % Plot steering drive for each condition
    plot(visObjPosition, steering_drive(1, :), 'DisplayName', comparisonLabel{1}, 'LineWidth', 1.5);
    plot(visObjPosition, steering_drive(2, :), 'DisplayName', comparisonLabel{2}, 'LineWidth', 1.5);
    
    % Customize plot
    xlabel('Azimuthal Object Position (degrees)');
    ylabel('Steering Drive (DNa02R - DNa02L)');
    title(['Steering Drive vs. Object Position (' comparisonLabel{1} ' vs ' comparisonLabel{2} ')']);
    xlim([-180 180]);
    legend('show');
    grid on;
    hold off;
    
    % Save the plot as PNG and SVG
    saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_SteeringDrive.png']));
    set(gcf, 'renderer', 'Painters');  % Set renderer to 'Painters' for better vector graphic rendering
    saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_SteeringDrive.svg']));
    
end
