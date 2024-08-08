% m is the lifetime activity of some cell population
function m=ELU(m)
    m=rescale(m,-1,1);
    m=(m).*(m >= 0) + (exp(m) - 1).*(m < 0);
    m=rescale(m,0,1); % rescale to range from 0 to 1
end