% adjELU
% This function applies the adjusted Exponential Linear Unit (ELU) transformation 
% to the input array `m`, which represents the lifetime activity of a cell population. 
% The transformation is controlled by the `alpha` parameter, which compresses the 
% nonlinear portion, and the `shift` parameter, which specifies where the nonlinear 
% portion begins. The input is first rescaled to a range of (-1, 1), transformed, 
% and then rescaled to a range of (0, 1).
%
% INPUTS:
%   m     - Array of lifetime activity values for a cell population.
%   alpha  - Controls the compression of the nonlinear portion.
%   shift  - Specifies where the nonlinear portion begins.
%
% OUTPUT:
%   m     - Transformed and rescaled array after applying the adjELU function.
%
function m = adjELU(m, alpha, shift)
    m = rescale(m, -1, 1);  % Rescale the input to a standard range (-1, 1)
    
    % Apply the ELU transformation with shift and alpha for compression
    m = (m - shift).*(m >= shift) + alpha * (exp(m - shift) - 1).*(m < shift);
    
    m = rescale(m, 0, 1);  % Rescale the output to range from 0 to 1
end
