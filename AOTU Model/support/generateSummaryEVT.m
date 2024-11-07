% generateSummaryEVT.m
%
% This function generates summary plots for average angular velocity versus object position
% across different steering gain (k) values, allowing a comparison between two conditions.
% Each k value is displayed in a separate tile, and both conditions are plotted together.
%
% INPUTS:
%   nK                - Number of steering gain values (k) in the experiment.
%   kValues           - Vector of k values used in the experiment.
%   evt_1             - Matrix of average angular velocity for Condition 1 across position bins.
%   evt_2             - Matrix of average angular velocity for Condition 2 across position bins.
%   posBins           - Vector defining the position bins (degrees).
%   comparisonLabel   - Cell array with labels for the two conditions being compared.
%   folder            - Folder directories
%
% OUTPUTS:
%   Generates a figure with plots showing the relationship between angular velocity and position for each k value.
%
% CREATED: 10/30/2024 - MC

function generateSummaryEVT(nK, kValues, evt_1, evt_2, posBins, comparisonLabel, folder)

    % Generate summary plot for Angular Velocity vs Object Position
    figure;
    set(gcf, 'Position', [100 100 1500 500]);  % Set figure size
    tiledlayout(1, nK, 'TileSpacing', 'compact');  % One row for each k value

    % Loop over k values for plotting
    for kIdx = 1:nK
        nexttile; hold on;
        
        % Plot the relationship for both conditions
        plot(posBins, evt_1(kIdx, :), 'DisplayName', comparisonLabel{1});
        plot(posBins, evt_2(kIdx, :), 'DisplayName', comparisonLabel{2});

        % Customize the plot
        xlabel('Position (deg)');
        ylabel('Avg. Angular Velocity (deg/s)');
        title(['k=' num2str(kValues(kIdx))]);
        if kIdx == 1
            legend('show', 'Location', 'northwest');
        end
        grid on;
        ylim([-600 600])
        xlim([-100 100])
    end
    sgtitle(['Binned Angular Velocity vs Object Position for ' comparisonLabel{1} ' vs ' comparisonLabel{2} ' Across k Values']);
    % Save EVT plot in both PNG and SVG formats
    saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_EVT.png']));
    set(gcf, 'renderer', 'Painters');  % Ensure high-resolution SVG output by setting renderer
    saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_EVT.svg']));

    % Apply zoomed x-axis limits to each tile
    for tileIdx = 1:nK
        nexttile(tileIdx);  % Select each tile in the current figure
        xlim([-50 50]);     % Set x-axis limits for zoomed view
    end

    % Save zoomed EVT plot in both PNG and SVG formats
    saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_EVT_Zoomed.png']));
    set(gcf, 'renderer', 'Painters');  % Set renderer for vector quality in zoomed version
    saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_EVT_Zoomed.svg']));

end
