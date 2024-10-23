% construct_2d_path
% This function generates and plots 2D paths for two sets of visual 
% object history data, allowing for comparison across different conditions. 
% The paths are constructed by unwrapping the position data and 
% creating a constantly increasing y position.
%
% INPUTS:
%   visobj_history1 - Matrix containing position data for condition 1.
%   visobj_history2 - Matrix containing position data for condition 2.
%   plotColors - Cell array of colors for plotting each condition.
%   percentage - Percentage of data points to plot.
%   kValues - Array of steering gains for labeling each condition.
%   comparisonLabel - Cell array of labels for the two conditions being compared.
%
% OUTPUTS:
%   path1 - Cell array of 2D paths for condition 1.
%   path2 - Cell array of 2D paths for condition 2.
% 
% This function ensures that both conditions have the same number of 
% rows and columns, calculates the number of data points to plot, 
% and generates a plot with labeled axes and legends.
%
function [path1, path2] = construct_2d_path(visobj_history1, visobj_history2, plotColors, percentage, kValues, comparisonLabel)
    % Get the number of conditions for each history (rows)
    numConditions1 = size(visobj_history1, 1);
    numConditions2 = size(visobj_history2, 1);
    
    % Ensure both histories have the same number of conditions (rows)
    if numConditions1 ~= numConditions2
        error('Both conditions must have the same number of conditions (rows).');
    end
    
    % Get the number of data points (columns)
    numDataPoints1 = size(visobj_history1, 2);
    numDataPoints2 = size(visobj_history2, 2);
    
    % Ensure both histories have the same number of data points (columns)
    if numDataPoints1 ~= numDataPoints2
        error('Both conditions must have the same number of data points (columns).');
    end
    
    % Calculate the number of data points to plot based on the percentage
    numPointsToPlot = round(percentage / 100 * numDataPoints1);
    
    % Generate a new figure and set its size
    figure;
    set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
    
    % Create a tiled layout for plotting
    tiledlayout(1, numConditions1);
    
    % Initialize cells to store paths for each condition
    path1 = cell(1, numConditions1);
    path2 = cell(1, numConditions1);
    
    % Loop through each condition (row)
    for cond = 1:numConditions1
        % Unwrap the position data between -180 and 180 for condition 1
        visobj1_unwrapped = unwrap(deg2rad(visobj_history1(cond, 1:numPointsToPlot)), pi);
        visobj1_unwrapped = rad2deg(visobj1_unwrapped);
        
        % Unwrap the position data between -180 and 180 for condition 2
        visobj2_unwrapped = unwrap(deg2rad(visobj_history2(cond, 1:numPointsToPlot)), pi);
        visobj2_unwrapped = rad2deg(visobj2_unwrapped);
        
        % Create a constantly increasing y position for both conditions
        y_positions = 1:numPointsToPlot;
        
        % Combine x and y positions into 2D paths
        path1{cond} = [visobj1_unwrapped(:), y_positions(:)];
        path2{cond} = [visobj2_unwrapped(:), y_positions(:)];
        
        % Plot both conditions in a comparison format
        nexttile;
        hold on;
        plot(path1{cond}(:, 1), path1{cond}(:, 2), 'Color', plotColors{1}, 'LineWidth', 2, 'DisplayName', comparisonLabel{1});
        plot(path2{cond}(:, 1), path2{cond}(:, 2), 'Color', plotColors{2}, 'LineWidth', 2, 'DisplayName', comparisonLabel{2});
        
        % Add labels, title, and legend
        xlabel('Unwrapped X Position (degrees)');
        ylabel('Y Position');
        title([num2str(kValues(cond)) ' X']);
        
        if cond == 1
            legend('show');
        end
        
        hold off;
    end
    
    % Ensure the x-axis is the same across all tiles
    linkaxes(findall(gcf, 'Type', 'axes'), 'x');
end
