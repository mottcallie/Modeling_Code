function [steerOutput] = dedicatedSystem(visualInput,rfSize,gain)

% remove values outside of receptive field
inputRestrict = visualInput;
inputRestrict(abs(visualInput)>rfSize) = 0;
x = inputRestrict;

steerOutput_raw = -(x+rfSize).*(x-rfSize).*(x);
steerOutput_adj = (steerOutput_raw./max(steerOutput_raw)).*rfSize;

steerOutput = steerOutput_adj*gain;
end

