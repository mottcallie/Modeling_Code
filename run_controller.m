% simple script to repeatedly run the controller simulation

% set controller properties
controllerDelay=200; %msec

% set folder
modelFolder = 'C:\Users\wilson\Dropbox (HMS)\MATLAB\Modeling';

%% simple sawtooth
cd(modelFolder)

% generate sawtooth stimulus
sweepRange = 75;    %degrees
sweepSpeed = 100;    %degrees/sec
duration = 10;      %seconds
funcFreq = 1000;

stimwave = generate_sawtooth(sweepRange,sweepSpeed,duration,funcFreq);
cc_pd = controlsimulation(stimwave,funcFreq,controllerDelay,1,1);
sgtitle(['Model Performance: target moving at ' num2str(sweepSpeed) 'mm/s'])

% save plot of trial data
cd([modelFolder '\plots'])
saveas(gcf,['model_' num2str(sweepRange) 'sweep_' num2str(sweepSpeed) 'degsec_' num2str(duration) 'sec' '_plot.png']);

%% pseudorandom target
cd(modelFolder)

% generate coherent path
sweepRange = 50;   %degrees
duration = 60;     %seconds
funcFreq = 1000;

stimwave = generate_coherentpath(duration,sweepRange,funcFreq);
cc_pd = controlsimulation(stimwave,funcFreq,controllerDelay,2,1);
sgtitle(['Model Performance: coherent path'])

% save plot of trial data
cd([modelFolder '\plots'])
saveas(gcf,['model_' 'coherentpath_plot.png']);