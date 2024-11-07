% ELU
% This function applies the Exponential Linear Unit (ELU) transformation 
% to a given input array `m`, which represents the lifetime activity of 
% a cell population. The function rescales the input to the range of 
% -1 to 1, applies the ELU transformation, and then rescales the result 
% to the range of 0 to 1.
%
% INPUTS:
%   m - Array of lifetime activity values for a cell population.
%
% OUTPUT:
%   m - Transformed and rescaled array after applying the ELU function.
%
function m=ELU(m)
    m=rescale(m,-1,1);
    m=(m).*(m >= 0) + (exp(m) - 1).*(m < 0);
    m=rescale(m,0,1); % rescale to range from 0 to 1
end