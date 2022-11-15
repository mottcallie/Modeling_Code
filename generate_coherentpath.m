function  stimwave = generate_coherentpath(duration,sweepRange,funcFreq)
%% important variables
% adjust variables to adjust path smoothing
smoothNoise = .1; %increase to increase gaussian smoothing

%% generate random white noise
b = 10;
buff = zeros(1,b);
nb = duration+b*2; %add buffer to start/end

% random integers within specified interval
%whitenoise = randi([amp*-1 amp],1,nb); %opt1: rand int
whitenoise = randn(1,duration)*sweepRange; %opt2: rand normal distribution
whitenoise = [buff whitenoise buff];

%% smooth white noise
smoothnoise=smoothdata(whitenoise,'gaussian',smoothNoise);

%% interpolate smoothed noise
x = 1:1:nb;
xi = 1/funcFreq:1/funcFreq:nb;

intnoise = interp1(x,smoothnoise,xi,'spline');
% interpolating often generates artifacts at start/end
%intnoise(1:b*funcFreq) = []; %remove buffer at start
%intnoise(end-b*funcFreq+1:end) = []; %remove buffer at end

%% store final path
stimwave = intnoise;

end

