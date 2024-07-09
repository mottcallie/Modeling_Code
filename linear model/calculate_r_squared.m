% calculate r2 value
% Y - actual response array
% Yfit - fit based on predictor
%
function r_squared = calculate_r_squared(observed, predicted)
    % Check if the input vectors have the same length
    if length(observed) ~= length(predicted)
        error('Input vectors must have the same length');
    end
    
    % Calculate mean of observed values
    mean_observed = mean(observed,'omitnan');
    
    % Calculate total sum of squares (TSS)
    total_sum_squares = sum((observed - mean_observed).^2,'omitnan');
    
    % Calculate sum of squared errors (SSE)
    sum_squared_errors = sum((observed - predicted).^2,'omitnan');
    
    % Calculate R-squared
    r_squared = 1 - (sum_squared_errors / total_sum_squares);
end