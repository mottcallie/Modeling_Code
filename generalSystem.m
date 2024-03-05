function [steerOutput] = generalSystem(visualInput,rfSize,gain)

% remove values outside of receptive field
inputRestrict = visualInput;
inputRestrict(abs(visualInput)>rfSize) = 0;

steerOutput = inputRestrict*gain;
end

