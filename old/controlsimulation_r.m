%controlsimulation
% model behavior of a PD controller when presented with an oscillating
% target
%
% INPUT
% sweepSpeed        - speed of target (e.g., 50 70 100deg/s)
% controllerDelay   - controller delay (e.g. 300msec.)
% plotSim           - 1 to plot simulation, 0 to not
%
% OUTPUT
% cc_pd             - correlation coefficient for model performance
%
% CREATED - 09/12/2022 MC and RW
% UPDATED - 09/19/2022 MC removed weight variable, adjusted derivative amp
%

function [cc_pd] = controlsimulation_r(sweepSpeed,controllerDelay,plotSim)
%% initialize

% generate figure
if plotSim
    figure(1); clf;
    set(gcf,'Position',[100 100 1500 800])
    lw = 1;
end

% set basic parameters
sweepRange = 110;   %degrees, rough range for object position
duration = 200;     %seconds

smoothNoise = .1;   %increase to increase gaussian smoothing
funcFreq= 1000;


%% generate waveform

b = 50;
buff = zeros(1,b);
nb = duration+b*2; %add buffer to start/end

% random integers within specified interval
%whitenoise = randi([amp*-1 amp],1,nb); %opt1: rand int
whitenoise = randn(1,duration)*sweepRange; %opt2: rand normal distribution
whitenoise = [buff whitenoise buff];

% smooth white noise
smoothnoise=smoothdata(whitenoise,'gaussian',smoothNoise);

% interpolate smoothed noise
x = 1:1:nb;
xi = 1/funcFreq:1/funcFreq:nb;

intnoise = interp1(x,smoothnoise,xi,'spline');
% interpolating often generates artifacts at start/end
intnoise(1:b*funcFreq) = []; %remove buffer at start
intnoise(end-b*funcFreq+1:end) = []; %remove buffer at end

% store final path
stimwave = intnoise;
t = 1/funcFreq:1/funcFreq:duration;

% plot target waveform
if plotSim
    s1 = subplot(2,1,1);
    plot(t,stimwave,'Color','k','LineWidth',lw)
    yline(0,'Color','#7E2F8E')
    ylabel('Target Position (deg)');
    axis tight
end


%% calculate error and output for each controller unit SEPERATELY

% generate delay array
delay = zeros(1,round((controllerDelay/1000)*funcFreq));

% generate position unit
pEC = stimwave;                      %position error = obj position
% add delay
pEC = [delay pEC];
pEC = pEC(1:length(t));
% normalize derivative error unit
%pEC = 2*amp*pEC/range(pEC);

% generate derivative error unit
dEC = diff([stimwave nan])/(1/funcFreq);  %derivative error = derivative/change t
% add delay
dEC = [delay dEC];
dEC = dEC(1:length(t));
% normalize derivative error unit
%dEC = 2*amp*dEC/range(dEC(3*sfq:7*sfq));

% plot position and derivative controller outputs
if plotSim
    s2 = subplot(2,1,2);
    hold on
    plot(t,pEC,'Color','#D95319','LineWidth',lw)
    plot(t,dEC,'Color','#77AC30','LineWidth',lw)
end


%% calculate COMBINED error controller

% sum PD control
pdEC = pEC + dEC;

% calculate correlation coefficient between waveform input and PD output
cc_pd = corrcoef(stimwave,pdEC);

% plot combined PD controller output
if plotSim
    plot(t,pdEC,'Color','#0072BD','LineWidth',lw*2)

    axis tight
    ylabel('Controller Output (deg/s)');
    %ylim([-sweepSpeed sweepSpeed])
    linkaxes([s1,s2],'x');
    
    % add labels
    sgtitle(['Model Performance: target moving at' num2str(sweepSpeed) 'mm/s'])
    xlabel('Time (s)')
end

end
