% generateSummaryInputDistributions
%
% This function generates summary plots comparing the distributions of input history 
% for two different conditions across multiple steering gain (k) values. Each k value 
% is displayed in a separate tile, allowing a direct comparison between the two conditions.
% Additionally, an input-output curve for DNa02 is plotted above each distribution comparison.
%
% INPUTS:
%   nK                - Number of steering gain values (k) in the experiment.
%   kValues           - Vector of k values used in the experiment.
%   inputDist_1       - Matrix of input history distributions for Condition 1 across bins.
%   inputDist_2       - Matrix of input history distributions for Condition 2 across bins.
%   bins              - Vector defining the bins for input values.
%   comparisonLabel   - Cell array with labels for the two conditions being compared.
%   runSettings       - Struct containing the parameters, including DNa02input and shift.
%   folder            - Folder directories for saving figures.
%
% OUTPUTS:
%   Generates a figure with plots showing the DNa02 input-output curves and 
%   distributions for each k value.
%
% CREATED: 11/05/2024 - MC

function generateSummaryInputDistributions(nK, kValues, inputDist_1, inputDist_2, bins, comparisonLabel, runSettings, folder)

    % Motor input-output for DNa02
    DNa02input = runSettings.DNa02input;  % Input range for downstream neurons
    alpha = 1; % Specify alpha based on model requirements
    DNa02output = adjELU(DNa02input, alpha, runSettings.shift); % Output range for downstream neurons

    % Generate summary plot for DNa02 Input-Output Curve and Input History Distributions
    figure;
    set(gcf, 'Position', [100 100 1400 300]);  % Set figure size
    tiledlayout(2, nK, 'TileSpacing', 'compact');  % Two rows for each k value
    xrange = [runSettings.minIn runSettings.maxIn];
    % Plot DNa02 input-output curve in the first row
    for kIdx = 1:nK
        nexttile(kIdx); hold on;
        
        % Plot DNa02 input-output curve
        plot(DNa02input, DNa02output, 'DisplayName', 'DNa02 Input-Output Curve');

        % Customize the plot
        xlabel('Input Value');
        ylabel('Output Value');
        title(['k=' num2str(kValues(kIdx))]);
        xline(0,'k')
        axis tight
        xlim(xrange)
        grid on;
    end

    % Plot input history distributions for both conditions in the second row
    for kIdx = 1:nK
        nexttile(kIdx + nK); hold on;
        
        % Plot the distributions for both conditions
        plot(bins, inputDist_1(kIdx, :), 'DisplayName', comparisonLabel{1});
        plot(bins, inputDist_2(kIdx, :), 'DisplayName', comparisonLabel{2});

        % Customize the plot
        xlabel('Input Value Bins');
        ylabel('Normalized Frequency');
        if kIdx == 1
            legend('show', 'Location', 'northwest');
        end
        grid on;
        xline(0,'k')
        xlim(xrange)
        ylim([0 1]);  % Set y-axis to normalized range
    end
    sgtitle(['DNa02 Input-Output and Input History Distribution Comparison: ' comparisonLabel{1} ' vs ' comparisonLabel{2} ' Across k Values']);

    % Save the plot in both PNG and SVG formats
    saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_Input_Distributions_with_DNa02_Curve.png']));
    set(gcf, 'renderer', 'Painters');  % Ensure high-resolution SVG output by setting renderer
    saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_Input_Distributions_with_DNa02_Curve.svg']));

end
