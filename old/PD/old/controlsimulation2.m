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
% bestWeight        - weight of derivative error controller that results in
%                     best performance of the PD output
% bestPerfm         - correlation coefficient for best model performance
%
% CREATED - 09/12/2022 MC and RW
% UPDATED - 09/16/2022 MC model selects D weight based on correlation
%

function [bestWeight,bestPerfm] = controlsimulation(sweepSpeed,controllerDelay,plotSim)
%% initialize

% generate figure
if plotSim
    figure(1); clf;
    set(gcf,'Position',[100 100 1500 800])
    lw = 1;
end

% set basic parameters
sweepRange = 70;    %degrees
duration = 10;      %seconds
sfq = 1000;         %sample frequency


%% generate waveform
% note: also possible with sawtooth function, but it is annoying to use

% create single left/right sweep
amp = sweepRange/2;
sweepAway = 0:1:amp;
sweepBack = amp-1:-1:1;
sweep = [(0+sweepAway) (0+sweepBack) (0-sweepAway) (0-sweepBack)]; %combine

% stretch sweep according to sweep speed
speed_sfq = round(sfq/sweepSpeed);
stimwave_single = repelem(sweep,speed_sfq);
stimwave_dur = length(stimwave_single);

% create start/stop buffer
buffer = zeros(1,(1 * sfq));    %seconds * sample frequency

% generate full target waveform
duration_sfq = duration*sfq;
nsweeps = round(duration_sfq/stimwave_dur);         %rough # of sweeps given duration
stimwave_multi = repmat(stimwave_single,1,nsweeps); %repeat sweeps
stimwave = [buffer stimwave_multi buffer];          %add in buffers

% set time
d = length(stimwave)/sfq; %duration of stimwave (sec.)
t = 1/sfq:1/sfq:d;

% plot target waveform
if plotSim
    s1 = subplot(2,1,1);
    plot(t,stimwave,'Color','k','LineWidth',lw)
    yline(0,'Color','#7E2F8E')
    ylabel('Target Position (deg)');
    axis tight
end


%% calculate error and output for each controller unit SEPERATELY

% smooth stimwave input
sm_stimwave = smoothdata(stimwave,'gaussian',1.2*sfq);
% generate delay array
delay = zeros(1,(controllerDelay/1000)*sfq);

% generate position unit
pEC = sm_stimwave;              %sposition error = obj position
% add delay
pEC = [delay pEC];
pEC = pEC(1:length(t));
% normalize derivative error unit
pEC = 2*amp*pEC/range(pEC);

% generate derivative error unit
dEC = diff([sm_stimwave nan]);  %derivative error = first derivative
% add delay
dEC = [delay dEC];
dEC = dEC(1:length(t));
% normalize derivative error unit
dEC = 2*amp*dEC/range(dEC(3*sfq:7*sfq));

% plot position and derivative controller outputs
if plotSim
    s2 = subplot(2,1,2);
    hold on
    plot(t,pEC,'Color','#D95319','LineWidth',lw)
    plot(t,dEC,'Color','#77AC30','LineWidth',lw)
end


%% calculate COMBINED error controller

% set weight of P vs D controller units based on model performance
% Note: in this case, we are basing "performance" on input/output correlation
weightRange = 0:0.1:2;      %range of weights to test
weightPerfm = weightRange'; %load in for later
for nd = 1:length(weightRange)
    % select D weight
    tryWeight = weightRange(nd);
    % sum PD control
    pdEC(nd,:) = pEC + (tryWeight*dEC);
    % normalize PD controller
    pdEC(nd,:) = 2*amp*pdEC(nd,:)/range(pdEC(nd,:));

    % calculate correlation coefficient between waveform input and PD output
    cc_pd = corrcoef(stimwave,pdEC(nd,:));
    weightPerfm(nd,2) = cc_pd(2);
end

% evaluate best weight, based on best model performance
[bestPerfm,r] = max(weightPerfm(:,2));
bestWeight = weightPerfm(r,1);

% select best PD controller based on D weight
pdECf = pdEC(r,:);

% plot combined PD controller output
if plotSim
    plot(t,pdECf,'Color','#0072BD','LineWidth',lw*2)

    axis tight
    ylabel('Controller Output (deg/s)');
    xlabel('Time (s)')
    ylim([-amp amp])
    linkaxes([s1,s2],'x');
end

end
