function stimwave = generate_sawtooth(sweepRange,sweepSpeed,duration,funcFreq)
% create single left/right sweep
amp = sweepRange/2;
sweepAway = 0:1:amp;
sweepBack = flip(sweepAway(2:end-1));
sweep = [(0+sweepAway) (0+sweepBack) (0-sweepAway) (0-sweepBack)]; %combine

% stretch sweep according to sweep speed
speed_sfq = round(funcFreq/sweepSpeed);
stimwave_single = repelem(sweep,speed_sfq);
stimwave_dur = length(stimwave_single);

% create start/stop buffer
buffer = zeros(1,(1 * funcFreq));    %seconds * sample frequency

% generate full target waveform
duration_sfq = duration*funcFreq;
nsweeps = round(duration_sfq/stimwave_dur);         %rough # of sweeps given duration
stimwave_multi = repmat(stimwave_single,1,nsweeps); %repeat sweeps
stimwave = [buffer stimwave_multi buffer];          %add in buffers
end