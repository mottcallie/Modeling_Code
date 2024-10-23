function indices = findFirstZeroCrossingMultiple(traces)

    % Initialize the output vector
    indices = zeros(size(traces, 1), 1);
    maxIdx = size(traces,2);

    % Loop through each trace
    for i = 1:size(traces, 1)
        trace = traces(i, :);
        % Find indices where the sign changes
        signChanges = find(diff(sign(trace)));
        % Determine the first zero crossing
        if ~isempty(signChanges)
            indices(i) = signChanges(1) + 1;
        else
            indices(i) = maxIdx; % Return max if there are no zero crossings
        end
    end
end