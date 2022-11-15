%controlsimulation
% model behavior of a PD controller when presented with an oscillating
% target
%
% INPUT
% oSpeed            - speed of target (e.g., 50 70 100deg/s)
% controllerDelay   - controller delay (e.g. 300msec.)
%
% CREATED - MC and RW 09/12/2022
%

function controlsimulation(oSpeed,controllerDelay)
%% initialize
% generate figure
figure(1); clf;
set(gcf,'Position',[100 100 1500 800])
lw = 1;


%% generate waveform
% set waveform variables
simDuration = 10;   %duration (sec.)
amp = 35;           %amplitude
ffq = 1;            %fundamental frequency
fs = 1000;          %sample frequency

y = oSpeed/(amp*2)/2; %stimulus velocity (convert to Hz)
T = simDuration*(1/ffq);

% generate target waveform
t = 0:1/fs:T-1/fs;  %time series
stimwave = sawtooth(2*pi*y*t,1/2);
wavestart = find(stimwave==0,1);                %find first midline cross
stimwave = circshift(stimwave,-wavestart+1,2);  %shift to start at 0
stimwave = stimwave*amp;                        %multiply waveform by amp
% add buffer to target waveform
nb = 1*fs;                                       %buffer length (msec)
stimwave = [zeros(1,nb) stimwave zeros(1,nb)];  %add start/stop buffer
t = [t (t(1:nb*2)+t(end))];                     %extend time series

% plot target waveform
s(1) = subplot(2,1,1);
plot(t,stimwave,Color='k',LineWidth=lw)
yline(0,'Color','#7E2F8E')
s(1) = ylabel('Target Position (deg)');
axis tight


%% calculate error and output for each controller unit SEPERATELY

% smooth stimwave input
sm_stimwave = smoothdata(stimwave,'gaussian',1.2*fs);
% generate delay array
delay = zeros(1,(controllerDelay/1000)*fs);

% generate position unit
pEC = sm_stimwave;              %sposition error = obj position
% add delay
pEC = [delay pEC];
pEC = pEC(1:length(t));
% normalize derivative error unit
pEC = 2*amp*pEC/range(pEC(3*fs:7*fs));

% generate derivative error unit
dEC = diff([sm_stimwave nan]);  %derivative error = first derivative
% add delay
dEC = [delay dEC];
dEC = dEC(1:length(t));
% normalize derivative error unit
dEC = 2*amp*dEC/range(dEC(3*fs:7*fs));

% plot position and derivative controller outputs
s(2) = subplot(2,1,2);
hold on
plot(t,pEC,Color='#D95319',LineWidth=lw)
plot(t,dEC,Color='#77AC30',LineWidth=lw)


%% calculate COMBINED error controller

% set how much to weight derivative error unit
% >1 will favor the derivative error unit, <1 will favor the position unit
dWeight = 1+y;

% sum PD control
pdEC = pEC + (dWeight*dEC);
% normalize PD controller
pdEC = 2*amp*pdEC/range(pdEC(3*fs:7*fs));

% plot combined PD controller output
plot(t,pdEC,Color='#0072BD',LineWidth=lw*2)

axis tight
xlabel('Time (s)')
ylim([-amp amp])
linkaxes(s,'x');


end
