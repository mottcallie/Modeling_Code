% plot_adjELU
% This function plots the adjusted Exponential Linear Unit (adjELU) transformation 
% for different alpha values.
%
% INPUTS:
%   alpha_values - A vector of alpha values to compare.
%   runSettings - Structure containing input values, specifically in runSettings.DNa02input,
%                 and the shift value for the adjELU function.
%
% OUTPUTS:
%   Generates a plot of the transformed values for each alpha value against the input values.
%
function plot_adjELU(alpha_values, runSettings)
    % Use input values from runSettings.DNa02input
    input_values = runSettings.DNa02input;
    shift = runSettings.shift;

    % Initialize figure
    hold on;

    % Loop through different alpha values and plot results
    for i = 1:length(alpha_values)
        alpha = alpha_values(i);
        output_values = adjELU(input_values, alpha, shift);
        
        % Plot the transformed values for each alpha
        plot(input_values, output_values, 'DisplayName', ['\alpha = ', num2str(alpha)],'LineWidth',1);
    end

    % Customize plot
    xlabel('Input values');
    ylabel('Transformed values');
    title('ELU');
    legend('show','Location','northwest');
    grid on;
    axis tight
    hold off;
    
end
