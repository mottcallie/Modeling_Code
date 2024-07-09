% simple script to repeatedly run the controller simulation

% set controller properties


% set folder
modelFolder = 'C:\Users\wilson\Dropbox (HMS)\MATLAB\Modeling_Code';

%% controller simulation in which target speed is varied
cd([modelFolder '\plots'])

% set stimulus parameters
sweepSpeeds = 50:20:170; %deg/sec
sweepRange = 75; %deg
duration = 10; %sec
funcFreq = 1000;
ctrlDelay = 200; %msec

% initialize
pAmp = [];
pPhase = [];
dAmp = [];
dPhase = [];

% for each sweep speed
for swp = 1:length(sweepSpeeds)
    thisSpeed = sweepSpeeds(swp);
    stimwave = generate_sawtooth(sweepRange,thisSpeed,duration,funcFreq);
    % run the controller simulation
    [~,amp,phase] = controlsimulation(stimwave,funcFreq,ctrlDelay,1,1);
    sgtitle(['target speed = ' num2str(thisSpeed) 'deg/s'])

    % store variables
    pAmp(swp) = amp.pos;
    pPhase(swp) = phase.pos;
    dAmp(swp) = amp.der;
    dPhase(swp) = phase.der;
    cAmp(swp) = amp.combo;
    cPhase(swp) = phase.combo;

    % save plot of trial data
    saveas(gcf,['model_' sprintf('%03d',thisSpeed) 'degsec' '_plot.png']);


    % optional wait to view each plot
    %wait = input('Next? ');
end


%% plot results
% plot results of position only controller
% initialize
cd([modelFolder '\plots'])
figure(2);clf;
set(gcf,'Position',[100 100 500 400])

% plot
subplot(1,2,1)
f1 = fit(sweepSpeeds',pAmp','poly2'); %fit
plot(f1,sweepSpeeds,pAmp,'o')
ylabel('amp (deg/s)')
xlabel('target speed (deg/sec)')
legend('off')
axis tight
ylim([0 80])

subplot(1,2,2)
f2 = fit(sweepSpeeds',pPhase','poly1'); %fit
plot(f2,sweepSpeeds,pPhase,'o')
ylabel('lag (msec)')
xlabel('target speed (deg/sec)')
legend('off')
axis tight
ylim([0 300])

% save
sgtitle('position controller only')
saveas(gcf,'summary_positioncontroller_speedbattery_plot.png');


% plot results of derivative only controller
% initialize
cd([modelFolder '\plots'])
figure(3);clf;
set(gcf,'Position',[100 100 500 400])

% plot
subplot(1,2,1)
f1 = fit(sweepSpeeds',dAmp','poly2'); %fit
plot(f1,sweepSpeeds,dAmp,'o')
ylabel('amp (deg/s)')
xlabel('target speed (deg/sec)')
legend('off')
axis tight
ylim([0 80])

subplot(1,2,2)
f2 = fit(sweepSpeeds',dPhase','exp1'); %fit
plot(f2,sweepSpeeds,dPhase,'o')
ylabel('lag (msec)')
xlabel('target speed (deg/sec)')
legend('off')
axis tight
ylim([-600 0])

% save
sgtitle('derivative controller only')
saveas(gcf,'summary_derivativecontroller_speedbattery_plot.png');


% plot results of position+derivative combined controller
% initialize
cd([modelFolder '\plots'])
figure(4);clf;
set(gcf,'Position',[100 100 500 400])

% plot
subplot(1,2,1)
f1 = fit(sweepSpeeds',cAmp','poly2'); %fit
plot(f1,sweepSpeeds,cAmp,'o')
ylabel('amp (deg/s)')
xlabel('target speed (deg/sec)')
legend('off')
axis tight
ylim([0 80])
subplot(1,2,2)
f2 = fit(sweepSpeeds',cPhase','poly2'); %fit
plot(f2,sweepSpeeds,cPhase,'o')
ylabel('lag (msec)')
xlabel('target speed (deg/sec)')
legend('off')
axis tight
ylim([-150 150])
yline(0)
% save
sgtitle('combined P+D controller')
saveas(gcf,'summary_combinedcontroller_speedbattery_plot.png');



