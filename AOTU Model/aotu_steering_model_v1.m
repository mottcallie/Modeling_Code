%% initialize
clear
close all

% pull file path
thisFile = matlab.desktop.editor.getActiveFilename;
[filePath,~,~] = fileparts(thisFile);
% jump to that folder
cd([filePath '/plots'])

visobj_position = linspace(-180, 180, 361); % visual object positions in azimuthal space ranges from -180 deg to +180 deg

%% generate neuron models

% generate visual tuning curve for AOTU019
curve_width = "thin";
switch curve_width
    case "wide"
        % visual tuning curves of each cell, over 360 deg azimuthal space
        AOTU019R=cos(deg2rad(visobj_position-35))';
        AOTU019L=cos(deg2rad(visobj_position+35))';
    case "thin"
        % visual tuning curves of each cell, over 360 deg azimuthal space
        f19 = 3; %frequency modifier
        AOTU019R=cos(f19*deg2rad(visobj_position-35))';
        AOTU019L=cos(f19*deg2rad(visobj_position+35))';
end
AOTU019R = ELU(AOTU019R); % output range is now [0 1]
AOTU019L = ELU(AOTU019L);
% restrict tuning curve to front field of view
if curve_width == "thin"
    iR = find(AOTU019R == 0);
    AOTU019R(1:iR(f19-1)) = 0;
    AOTU019R(iR(f19):end) = 0;

    iL = find(AOTU019L == 0);
    AOTU019L(1:iL(1)) = 0;
    AOTU019L(iL(f19-1):end) = 0;
end

% generate visual tuning curve for AOTU025
switch curve_width
    case "wide"
        AOTU025R=cos(deg2rad(visobj_position-150))';
        AOTU025L=cos(deg2rad(visobj_position+150))';
    case "thin"
        f25 = 2; %frequency modifier
        AOTU025R=cos(f25*deg2rad(visobj_position-100))';
        AOTU025L=cos(f25*deg2rad(visobj_position+100))';
end
AOTU025R = ELU(AOTU025R);
AOTU025L = ELU(AOTU025L);
% restrict tuning curve to peripheral field of view
if curve_width == "thin"
    iR = find(AOTU025R == 0);
    AOTU025R(1:iR(f25)) = 0;

    iL = find(AOTU025L == 0);
    AOTU025L(iL(f25-1):end) = 0;
end

% generate DNa02 output curve
DNa02input = linspace(-1, 2, 1000); % range of possible input values for DNa02
DNa02output = ELU(DNa02input); % range of possible output values for DNa02

% plot
figure
set(gcf,'Position',[100 100 500 800]);
tiledlayout(3,1,'TileSpacing','compact')
nexttile
plot(visobj_position,AOTU019R)
hold on
plot(visobj_position,AOTU019L)
ylabel("firing rate")
legend("AOTU019R","AOTU019L")
box off,xlim([-180 180]),xticks([-180 0 180])
nexttile
plot(visobj_position,AOTU025R)
hold on
plot(visobj_position,AOTU025L)
ylabel("firing rate")
legend("AOTU025R","AOTU025L")
box off,xlim([-180 180]),xticks([-180 0 180])
xlabel("visual object position (deg)")
nexttile
plot(DNa02input,DNa02output)
box off, xticks([-1 0 1 2]), yticks([0 1])
xlabel("DNa02 input"),ylabel("DNa02 output")
% save
saveas(gcf,'aotu_steering_model_inputs.png');


%% simulate pursuit #1

% set run options (e.g., inhibit, excite, slow)
AOTU019_manip = {"inhibitory"; "excitatory"; "slow"};
numCond = length(AOTU019_manip);

numTimepoints=60;
k=180; % free parameter specifying the maximum change in head direction per time step
noiseLevel = 10;

fs = 10; % simulation update rate, in Hz
fpass = 2; % lowpass cutoff for random component, in Hz
timebase=linspace(0,numTimepoints/fs,numTimepoints);

numRuns=50; % number of times to run the simulation
numPlots = 5; % number of runs to plot

random=zeros(numRuns,numTimepoints); % random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
startSide = ones(numRuns,1); %alternate which side the target starts on
startSide(2:2:end) = -1;

DNa02R_history=zeros(numRuns,numTimepoints); % DNa02R activity over time
DNa02L_history=zeros(numRuns,numTimepoints); % DNa02L activity over time
rot_vel_history=zeros(numRuns,numTimepoints); % fly's rotational velocity over time
visobj_history=zeros(numRuns,numTimepoints); % position of the visual object over time

% initialize
figure
set(gcf,'Position',[100 100 700 800]);
tiledlayout(numCond,1,'TileSpacing','compact')
all_runs = [];

% for each condition
for nc = 1:numCond
    % generated repeated runs
    for r=1:numRuns
        rng(r+2) % initialize random number generator to produce a new frozen noise snippet
        random(r,:) = randn(numTimepoints,1); % drawn from a Gaussian
        random(r,:) = lowpass(random(r,:),fpass,fs); % low-pass filter
        random(r,:) = noiseLevel*zscore(random(r,:)); % rescale
        visobj_history(r,1:2)=130*startSide(r); % starting visual object position, in degrees
        for t=3:numTimepoints
            v_tminus1=find(visobj_position==wrapTo180(round(visobj_history(r,t-1)))); % get the index value (v) of the visual object position at time t-1
            v_tminus2=find(visobj_position==wrapTo180(round(visobj_history(r,t-2)))); % get the index value (v) of the visual object position at time t-2

            switch AOTU019_manip{nc}
                case "inhibitory"
                    current_inputR = AOTU025R(v_tminus2) - AOTU019L(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) - AOTU019R(v_tminus1);
                case "excitatory"
                    current_inputR = AOTU025R(v_tminus2) + AOTU019R(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) + AOTU019L(v_tminus1);
                case "slow"
                    current_inputR = AOTU025R(v_tminus2) - AOTU019L(v_tminus2);
                    current_inputL = AOTU025L(v_tminus2) - AOTU019R(v_tminus2);
            end

            d = interp1(DNa02input,1:length(DNa02input),current_inputR,'nearest');
            DNa02R_history(r,t) = DNa02output(d);
            d = interp1(DNa02input,1:length(DNa02input),current_inputL,'nearest');
            DNa02L_history(r,t) = DNa02output(d);
            rot_vel_history(r,t) = k*(DNa02R_history(r,t) - DNa02L_history(r,t)) + random(r,t); % update steering command;
            visobj_history(r,t)=visobj_history(r,t-1) - rot_vel_history(r,t); % update visual object position (so that positive rotational velocity moves the object leftward)
            visobj_history(r,t)=wrapTo180(visobj_history(r,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
        end

    end

    % generate plot
    nexttile
    plot(timebase(1,2:end),visobj_history(1:numPlots,2:end),'--o');
    ylim([-180 180]), yticks([-180 0 180]), xticks([])
    yline(0), box off
    title(AOTU019_manip{nc})
    xlabel('time')
    ylabel('target pos (degrees)')

    % store results for additional analysis
    all_runs(:,:,nc) = visobj_history;
end
% save
saveas(gcf,'aotu_steering_model_simulationOutput.png');

% plot heading distribution
figure
set(gcf,'Position',[100 100 600 800]);
tiledlayout(numCond,2,'TileSpacing','compact')
a = 0.4; %opacity
for nc = 1:numCond
    nexttile
    % polar plot
    polarhistogram(deg2rad(all_runs(:,:,nc)),20,'FaceColor','#77AC30','FaceAlpha',a)
    title(AOTU019_manip{nc})

    nexttile
    %histogram
    histogram(all_runs(:,:,nc),'BinWidth',5,'FaceColor','#77AC30','FaceAlpha',a)
    title(AOTU019_manip{nc})
    xlabel('target pos (deg)')
    ylabel('counts')
end
% save
saveas(gcf,'aotu_steering_model_hist.png');


%% generate modified AOTU019 neuron models

% set run options
AOTU019_manip = {"overlap"; "less overlap"; "no overlap"; "wide"};
numCond = length(AOTU019_manip);

% generate shifted visual tuning curves
AOTU019R_shift1 = circshift(AOTU019R,10);
AOTU019L_shift1 = circshift(AOTU019L,-10);
AOTU019R_shift2 = circshift(AOTU019R,20);
AOTU019L_shift2 = circshift(AOTU019L,-20);
AOTU019R_shift3 = circshift(AOTU019R,30);
AOTU019L_shift3 = circshift(AOTU019L,-30);

% plot
figure
set(gcf,'Position',[100 100 500 800]);
tiledlayout(numCond,1,'TileSpacing','compact')
nexttile
plot(visobj_position,AOTU019R)
hold on
plot(visobj_position,AOTU019L)
ylabel("firing rate")
legend("AOTU019R","AOTU019L")
title(AOTU019_manip{1})
box off,xlim([-180 180]),xticks([-180 0 180])
nexttile
plot(visobj_position,AOTU019R_shift1)
hold on
plot(visobj_position,AOTU019L_shift1)
ylabel("firing rate")
legend("AOTU019R","AOTU019L")
title(AOTU019_manip{2})
box off,xlim([-180 180]),xticks([-180 0 180])
nexttile
plot(visobj_position,AOTU019R_shift2)
hold on
plot(visobj_position,AOTU019L_shift2)
ylabel("firing rate")
legend("AOTU019R","AOTU019L")
title(AOTU019_manip{3})
box off,xlim([-180 180]),xticks([-180 0 180])
nexttile
plot(visobj_position,AOTU019R_shift3)
hold on
plot(visobj_position,AOTU019L_shift3)
ylabel("firing rate")
legend("AOTU019R","AOTU019L")
title(AOTU019_manip{4})
box off,xlim([-180 180]),xticks([-180 0 180])
% save
saveas(gcf,'o_aotu_steering_model_inputs.png');

%% simulate pursuit #2

numRuns=50; % number of times to run the simulation
numPlots = 5; % number of runs to plot

numTimepoints=60;
k=180; % free parameter specifying the maximum change in head direction per time step
noiseLevel = 10;

fs = 10; % simulation update rate, in Hz
fpass = 2; % lowpass cutoff for random component, in Hz
timebase=linspace(0,numTimepoints/fs,numTimepoints);

random=zeros(numRuns,numTimepoints); % random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
startSide = ones(numRuns,1); %alternate which side the target starts on
startSide(2:2:end) = -1;

DNa02R_history=zeros(numRuns,numTimepoints); % DNa02R activity over time
DNa02L_history=zeros(numRuns,numTimepoints); % DNa02L activity over time
rot_vel_history=zeros(numRuns,numTimepoints); % fly's rotational velocity over time
visobj_history=zeros(numRuns,numTimepoints); % position of the visual object over time

% initialize
figure
set(gcf,'Position',[100 100 700 800]);
tiledlayout(numCond,1,'TileSpacing','compact')
all_runs = [];

% for each condition
for nc = 1:numCond
    % generated repeated runs
    for r=1:numRuns
        rng(r+2) % initialize random number generator to produce a new frozen noise snippet
        random(r,:) = randn(numTimepoints,1); % drawn from a Gaussian
        random(r,:) = lowpass(random(r,:),fpass,fs); % low-pass filter
        random(r,:) = noiseLevel*zscore(random(r,:)); % rescale
        visobj_history(r,1:2)=130*startSide(r); % starting visual object position, in degrees
        for t=3:numTimepoints
            v_tminus1=find(visobj_position==wrapTo180(round(visobj_history(r,t-1)))); % get the index value (v) of the visual object position at time t-1
            v_tminus2=find(visobj_position==wrapTo180(round(visobj_history(r,t-2)))); % get the index value (v) of the visual object position at time t-2

            switch AOTU019_manip{nc}
                case "overlap"
                    current_inputR = AOTU025R(v_tminus2) - AOTU019L(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) - AOTU019R(v_tminus1);
                case "less overlap"
                    current_inputR = AOTU025R(v_tminus2) - AOTU019L_shift1(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) - AOTU019R_shift1(v_tminus1);
                case "no overlap"
                    current_inputR = AOTU025R(v_tminus2) - AOTU019L_shift2(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) - AOTU019R_shift2(v_tminus1);
                case "wide"
                    current_inputR = AOTU025R(v_tminus2) - AOTU019L_shift3(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) - AOTU019R_shift3(v_tminus1);
            end

            d = interp1(DNa02input,1:length(DNa02input),current_inputR,'nearest');
            DNa02R_history(r,t) = DNa02output(d);
            d = interp1(DNa02input,1:length(DNa02input),current_inputL,'nearest');
            DNa02L_history(r,t) = DNa02output(d);
            rot_vel_history(r,t) = k*(DNa02R_history(r,t) - DNa02L_history(r,t)) + random(r,t); % update steering command;
            visobj_history(r,t)=visobj_history(r,t-1) - rot_vel_history(r,t); % update visual object position (so that positive rotational velocity moves the object leftward)
            visobj_history(r,t)=wrapTo180(visobj_history(r,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
        end

    end

    % generate plot
    nexttile
    plot(timebase(1,2:end),visobj_history(1:numPlots,2:end),'--o');
    ylim([-180 180]), yticks([-180 0 180]), xticks([])
    yline(0), box off
    title(AOTU019_manip{nc})
    xlabel('time')
    ylabel('target pos (degrees)')

    % store results for additional analysis
    all_runs(:,:,nc) = visobj_history;
end
% save
saveas(gcf,'o_aotu_steering_model_simulationOutput.png');

% plot heading distribution
figure
set(gcf,'Position',[100 100 600 800]);
tiledlayout(numCond,2,'TileSpacing','compact')
a = 0.4; %opacity
for nc = 1:numCond
    nexttile
    % polar plot
    polarhistogram(deg2rad(all_runs(:,:,nc)),20,'FaceColor','#77AC30','FaceAlpha',a)
    title(AOTU019_manip{nc})

    nexttile
    %histogram
    histogram(all_runs(:,:,nc),'BinWidth',5,'FaceColor','#77AC30','FaceAlpha',a)
    title(AOTU019_manip{nc})
    xlabel('target pos (deg)')
    ylabel('counts')
end
% save
saveas(gcf,'o_aotu_steering_model_hist.png');


random=zeros(numRuns,numTimepoints); % random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
startSide = ones(numRuns,1); %alternate which side the target starts on
startSide(2:2:end) = -1;

DNa02R_history=zeros(numRuns,numTimepoints); % DNa02R activity over time
DNa02L_history=zeros(numRuns,numTimepoints); % DNa02L activity over time
rot_vel_history=zeros(numRuns,numTimepoints); % fly's rotational velocity over time
visobj_history=zeros(numRuns,numTimepoints); % position of the visual object over time

% initialize
figure
set(gcf,'Position',[100 100 700 800]);
tiledlayout(numCond,1,'TileSpacing','compact')
all_runs = [];

% for each condition
for nc = 1:numCond
    % generated repeated runs
    for r=1:numRuns
        rng(r+2) % initialize random number generator to produce a new frozen noise snippet
        random(r,:) = randn(numTimepoints,1); % drawn from a Gaussian
        random(r,:) = lowpass(random(r,:),fpass,fs); % low-pass filter
        random(r,:) = noiseLevel*zscore(random(r,:)); % rescale
        visobj_history(r,1:2)=130*startSide(r); % starting visual object position, in degrees
        for t=3:numTimepoints
            v_tminus1=find(visobj_position==wrapTo180(round(visobj_history(r,t-1)))); % get the index value (v) of the visual object position at time t-1
            v_tminus2=find(visobj_position==wrapTo180(round(visobj_history(r,t-2)))); % get the index value (v) of the visual object position at time t-2

            switch AOTU019_manip{nc}
                case "overlap"
                    current_inputR = AOTU025R(v_tminus2) + AOTU019L(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) + AOTU019R(v_tminus1);
                case "less overlap"
                    current_inputR = AOTU025R(v_tminus2) + AOTU019R_shift1(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) + AOTU019L_shift1(v_tminus1);
                case "no overlap"
                    current_inputR = AOTU025R(v_tminus2) + AOTU019R_shift2(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) + AOTU019L_shift2(v_tminus1);
                case "wide"
                    current_inputR = AOTU025R(v_tminus2) + AOTU019R_shift3(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) + AOTU019L_shift3(v_tminus1);
            end

            d = interp1(DNa02input,1:length(DNa02input),current_inputR,'nearest');
            DNa02R_history(r,t) = DNa02output(d);
            d = interp1(DNa02input,1:length(DNa02input),current_inputL,'nearest');
            DNa02L_history(r,t) = DNa02output(d);
            rot_vel_history(r,t) = k*(DNa02R_history(r,t) - DNa02L_history(r,t)) + random(r,t); % update steering command;
            visobj_history(r,t)=visobj_history(r,t-1) - rot_vel_history(r,t); % update visual object position (so that positive rotational velocity moves the object leftward)
            visobj_history(r,t)=wrapTo180(visobj_history(r,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
        end

    end

    % generate plot
    nexttile
    plot(timebase(1,2:end),visobj_history(1:numPlots,2:end),'--o');
    ylim([-180 180]), yticks([-180 0 180]), xticks([])
    yline(0), box off
    title(AOTU019_manip{nc})
    xlabel('time')
    ylabel('target pos (degrees)')

    % store results for additional analysis
    all_runs(:,:,nc) = visobj_history;
end
% save
saveas(gcf,'oe_aotu_steering_model_simulationOutput.png');

% plot heading distribution
figure
set(gcf,'Position',[100 100 600 800]);
tiledlayout(numCond,2,'TileSpacing','compact')
a = 0.4; %opacity
for nc = 1:numCond
    nexttile
    % polar plot
    polarhistogram(deg2rad(all_runs(:,:,nc)),20,'FaceColor','#77AC30','FaceAlpha',a)
    title(AOTU019_manip{nc})

    nexttile
    %histogram
    histogram(all_runs(:,:,nc),'BinWidth',5,'FaceColor','#77AC30','FaceAlpha',a)
    title(AOTU019_manip{nc})
    xlabel('target pos (deg)')
    ylabel('counts')
end
% save
saveas(gcf,'oe_aotu_steering_model_hist.png');


%% simulate pursuit #3

% compare stability of making small error control excitatory vs inhibitory
AOTU019R_in = AOTU019R;
AOTU019L_in = AOTU019L;
AOTU019R_ex = AOTU019R_shift1;
AOTU019L_ex = AOTU019L_shift1;

AOTU019_manip = {"inhibitory"; "excitatory";"inhibitory"; "excitatory"};
numCond = length(AOTU019_manip);

numRuns=50; % number of times to run the simulation
numPlots = 5; % number of runs to plot
noiseMin = 10;
noiseMax = 40;
noiseLevel = [noiseMin, noiseMin, noiseMax, noiseMax];

random=zeros(numRuns,numTimepoints); % random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
startSide = ones(numRuns,1); %alternate which side the target starts on
startSide(2:2:end) = -1;

DNa02R_history=zeros(numRuns,numTimepoints); % DNa02R activity over time
DNa02L_history=zeros(numRuns,numTimepoints); % DNa02L activity over time
rot_vel_history=zeros(numRuns,numTimepoints); % fly's rotational velocity over time
visobj_history=zeros(numRuns,numTimepoints); % position of the visual object over time

% initialize
figure
set(gcf,'Position',[100 100 700 800]);
tiledlayout(numCond,1,'TileSpacing','compact')
all_runs = [];

% for each condition
for nc = 1:numCond
    % generated repeated runs
    for r=1:numRuns
        rng(r+2) % initialize random number generator to produce a new frozen noise snippet
        random(r,:) = randn(numTimepoints,1); % drawn from a Gaussian
        random(r,:) = lowpass(random(r,:),fpass,fs); % low-pass filter
        random(r,:) = noiseLevel(nc)*zscore(random(r,:)); % rescale
        visobj_history(r,1:2)=130*startSide(r); % starting visual object position, in degrees
        for t=3:numTimepoints
            v_tminus1=find(visobj_position==wrapTo180(round(visobj_history(r,t-1)))); % get the index value (v) of the visual object position at time t-1
            v_tminus2=find(visobj_position==wrapTo180(round(visobj_history(r,t-2)))); % get the index value (v) of the visual object position at time t-2

            switch AOTU019_manip{nc}
                case "inhibitory"
                    current_inputR = AOTU025R(v_tminus2) - AOTU019L_in(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) - AOTU019R_in(v_tminus1);
                case "excitatory"
                    current_inputR = AOTU025R(v_tminus2) - AOTU019L_ex(v_tminus1);
                    current_inputL = AOTU025L(v_tminus2) - AOTU019R_ex(v_tminus1);
            end

            d = interp1(DNa02input,1:length(DNa02input),current_inputR,'nearest');
            DNa02R_history(r,t) = DNa02output(d);
            d = interp1(DNa02input,1:length(DNa02input),current_inputL,'nearest');
            DNa02L_history(r,t) = DNa02output(d);
            rot_vel_history(r,t) = k*(DNa02R_history(r,t) - DNa02L_history(r,t)) + random(r,t); % update steering command;
            visobj_history(r,t)=visobj_history(r,t-1) - rot_vel_history(r,t); % update visual object position (so that positive rotational velocity moves the object leftward)
            visobj_history(r,t)=wrapTo180(visobj_history(r,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
        end

    end

    % generate plot
    nexttile
    plot(timebase(1,2:end),visobj_history(1:numPlots,2:end),'--o');
    ylim([-180 180]), yticks([-180 0 180]), xticks([])
    yline(0), box off
    title(AOTU019_manip{nc})
    xlabel('time')
    ylabel('target pos (degrees)')

    % store results for additional analysis
    all_runs(:,:,nc) = visobj_history;
end
% save
saveas(gcf,'noise_aotu_steering_model_simulationOutput.png');

% plot heading distribution
figure
set(gcf,'Position',[100 100 600 800]);
tiledlayout(numCond,2,'TileSpacing','compact')
a = 0.4; %opacity
for nc = 1:numCond
    nexttile
    % polar plot
    polarhistogram(deg2rad(all_runs(:,:,nc)),20,'FaceColor','#77AC30','FaceAlpha',a)
    title(AOTU019_manip{nc})

    nexttile
    %histogram
    histogram(all_runs(:,:,nc),'BinWidth',5,'FaceColor','#77AC30','FaceAlpha',a)
    title(AOTU019_manip{nc})
    xlabel('target pos (deg)')
    ylabel('counts')
end
% save
saveas(gcf,'noise_aotu_steering_model_hist.png');

%% expansion of noise simulation

numRuns=5; % number of times to run the simulation
noiseLevel = noiseMin:10:noiseMax; % set range of noise values
numCond = length(noiseLevels);

random=zeros(numRuns,numTimepoints); % random component injected into rotational velocity signal (modeling the influence of other steering drives and/or noise)
startSide = ones(numRuns,1); %alternate which side the target starts on
startSide(2:2:end) = -1;

DNa02R_history=zeros(numRuns,numTimepoints); % DNa02R activity over time
DNa02L_history=zeros(numRuns,numTimepoints); % DNa02L activity over time
rot_vel_history=zeros(numRuns,numTimepoints); % fly's rotational velocity over time
visobj_history=zeros(numRuns,numTimepoints); % position of the visual object over time

% for both excitatory and inhibitory
for nt = 1:2
    % for each condition
    for nc = 1:numCond
        % generated repeated runs
        for r=1:numRuns
            rng(r+2) % initialize random number generator to produce a new frozen noise snippet
            random(r,:) = randn(numTimepoints,1); % drawn from a Gaussian
            random(r,:) = lowpass(random(r,:),fpass,fs); % low-pass filter
            random(r,:) = noiseLevel(nc)*zscore(random(r,:)); % rescale
            visobj_history(r,1:2)=130*startSide(r); % starting visual object position, in degrees
            for t=3:numTimepoints
                v_tminus1=find(visobj_position==wrapTo180(round(visobj_history(r,t-1)))); % get the index value (v) of the visual object position at time t-1
                v_tminus2=find(visobj_position==wrapTo180(round(visobj_history(r,t-2)))); % get the index value (v) of the visual object position at time t-2

                switch nt
                    case 1
                        current_inputR = AOTU025R(v_tminus2) - AOTU019L_in(v_tminus1);
                        current_inputL = AOTU025L(v_tminus2) - AOTU019R_in(v_tminus1);
                    case 2
                        current_inputR = AOTU025R(v_tminus2) - AOTU019L_ex(v_tminus1);
                        current_inputL = AOTU025L(v_tminus2) - AOTU019R_ex(v_tminus1);
                end

                d = interp1(DNa02input,1:length(DNa02input),current_inputR,'nearest');
                DNa02R_history(r,t) = DNa02output(d);
                d = interp1(DNa02input,1:length(DNa02input),current_inputL,'nearest');
                DNa02L_history(r,t) = DNa02output(d);
                rot_vel_history(r,t) = k*(DNa02R_history(r,t) - DNa02L_history(r,t)) + random(r,t); % update steering command;
                visobj_history(r,t)=visobj_history(r,t-1) - rot_vel_history(r,t); % update visual object position (so that positive rotational velocity moves the object leftward)
                visobj_history(r,t)=wrapTo180(visobj_history(r,t)); % ensure that the new visual object position lies within the range from -180 to +180 deg
            end

        end

        % pull distribution
        [counts,edges] = histcounts(visobj_history,'BinWidth',5)
    end
end

%% extra
% figure
% tiledlayout(2,1,'TileSpacing','compact')
% nexttile
% scatter(timebase(1,2:end),DNa02R_history(:,2:end),20)
% ylim([-1 2]), yticks([-1 0 1 2]), xticks([])
% yline(0), box off
% ylabel('DNa02R firing rate')
% nexttile
% scatter(timebase(1,2:end),DNa02L_history(:,2:end),20)
% ylim([-1 2]), yticks([-1 0 1 2]), xticks([])
% yline(0), box off
% ylabel('DNa02L firing rate')
% xlabel('time')

%% functions

function m=ELU(m) % m is the lifetime activity of some cell population
m=rescale(m,-1,1);
m=(m).*(m >= 0) + (exp(m) - 1).*(m < 0);
m=rescale(m,0,1); % rescale to range from 0 to 1
end

