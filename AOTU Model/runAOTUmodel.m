%% initialize
clear
close all

% pull file path
thisFile = matlab.desktop.editor.getActiveFilename;
[filePath,~,~] = fileparts(thisFile);
% jump to that folder
cd([filePath '/plots'])

% general simulation variables
visobj_position = linspace(-180, 180, 361); % visual object positions in azimuthal space ranges from -180 deg to +180 deg
nTime = 60; %time of simulation run
nPlot = 5; %number of example runs to plot

% plot variables
a = 0.4; %opacity
nPolarBins = 20; %polar histogram number of bins
histBinSize = 5; %normal histogram size of bins


%% vary AOTU019 synapse
% set variables
noiseLevel = 10;
AOTU019synapse = {"inhibitory","excitatory","slow"};
AOTU019overlap = "normal";
nCond = length(AOTU019synapse);

% initialize
figure
set(gcf,'Position',[100 100 1200 800]);
tiledlayout(nCond,6,'TileSpacing','compact')
all_runs = [];

% for each overlap condition
for nc = 1:nCond
    % select this condition
    thisSynapse = AOTU019synapse{nc};

    % run this simulation condition
    disp(join(["Simulating " thisSynapse " synapse condition..."],""))
    [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(noiseLevel,thisSynapse,AOTU019overlap,nTime);
    disp("Simulating complete.")
    % store this simulation condition
    all_runs(:,:,nc) = visobj_history;

    % generate condition plot
    % plot visual tuning
    nexttile
    plot(visobj_position,AOTU019,"Color","#0072BD")
    hold on
    plot(visobj_position,AOTU025,"Color","#7E2F8E")
    ylabel("firing rate")
    box off,xlim([-180 180]),xticks([-180 0 180])
    % plot example simulation runs
    nexttile([1 3])
    plot(timebase(1,2:end),visobj_history(1:nPlot,2:end),'--o');
    ylim([-180 180]), yticks([-180 0 180]), xticks([])
    yline(0), box off
    title(join([thisSynapse " synapse"],""))
    xlabel('time')
    ylabel('target pos (degrees)')
    % plot polar histogram
    nexttile
    polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a)
    %plot normal histogram
    nexttile
    histogram(visobj_history(:,3:end),'BinWidth',histBinSize,'FaceColor','#77AC30','FaceAlpha',a)
    xlabel('target pos (deg)')
    ylabel('counts')
    xline(0)
    %calculate standard deviation
    thisSTD = std(reshape(visobj_history(:,3:end),[],1));
    title(join(["std = " num2str(round(thisSTD))],""))
end

% save
saveas(gcf,'synapse_AOTUmodel.png');

%% vary binocular overlap

% set variables
noiseLevel = 10;
AOTU019synapse = "inhibitory";
AOTU019overlap = {"normal","medium","low","none"};
nCond = length(AOTU019overlap);

% initialize
figure
set(gcf,'Position',[100 100 1200 800]);
tiledlayout(nCond,6,'TileSpacing','compact')
all_runs = [];

% for each overlap condition
for nc = 1:nCond
    % select this condition
    thisSynapse = AOTU019overlap{nc};

    % run this simulation condition
    disp(join(["Simulating " thisSynapse " overlap condition..."],""))
    [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(noiseLevel,AOTU019synapse,thisSynapse,nTime);
    disp("Simulating complete.")
    % store this simulation condition
    all_runs(:,:,nc) = visobj_history;

    % generate condition plot
    % plot visual tuning
    nexttile
    plot(visobj_position,AOTU019,"Color","#0072BD")
    hold on
    plot(visobj_position,AOTU025,"Color","#7E2F8E")
    ylabel("firing rate")
    box off,xlim([-180 180]),xticks([-180 0 180])
    % plot example simulation runs
    nexttile([1 3])
    plot(timebase(1,2:end),visobj_history(1:nPlot,2:end),'--o');
    ylim([-180 180]), yticks([-180 0 180]), xticks([])
    yline(0), box off
    title(join([thisSynapse " binocular overlap"],""))
    xlabel('time')
    ylabel('target pos (degrees)')
    % plot polar histogram
    nexttile
    polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a)
    nexttile
    %plot normal histogram
    histogram(visobj_history(:,3:end),'BinWidth',histBinSize,'FaceColor','#77AC30','FaceAlpha',a)
    xlabel('target pos (deg)')
    ylabel('counts')
    xline(0)
    %calculate standard deviation
    thisSTD = std(reshape(visobj_history(:,3:end),[],1));
    title(join(["std = " num2str(round(thisSTD))],""))
end

% save
saveas(gcf,'overlap_AOTUmodel.png');


% vary binocular + excitatory overlap
% set variables
AOTU019synapse = "excitatory";

% initialize
figure
set(gcf,'Position',[100 100 1200 800]);
tiledlayout(nCond,6,'TileSpacing','compact')
all_runs = [];

% for each overlap condition
for nc = 1:nCond
    % select this condition
    thisSynapse = AOTU019overlap{nc};

    % run this simulation condition
    disp(join(["Simulating " thisSynapse " overlap condition w/ excitatory synapse..."],""))
    [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(noiseLevel,AOTU019synapse,thisSynapse,nTime);
    disp("Simulating complete.")
    % store this simulation condition
    all_runs(:,:,nc) = visobj_history;

    % generate condition plot
    % plot visual tuning
    nexttile
    plot(visobj_position,AOTU019,"Color","#0072BD")
    hold on
    plot(visobj_position,AOTU025,"Color","#7E2F8E")
    ylabel("firing rate")
    box off,xlim([-180 180]),xticks([-180 0 180])
    % plot example simulation runs
    nexttile([1 3])
    plot(timebase(1,2:end),visobj_history(1:nPlot,2:end),'--o');
    ylim([-180 180]), yticks([-180 0 180]), xticks([])
    yline(0), box off
    title(join([thisSynapse " binocular overlap"],""))
    xlabel('time')
    ylabel('target pos (degrees)')
    % plot polar histogram
    nexttile
    polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a)
    nexttile
    %plot normal histogram
    histogram(visobj_history(:,3:end),'BinWidth',histBinSize,'FaceColor','#77AC30','FaceAlpha',a)
    xlabel('target pos (deg)')
    ylabel('counts')
    xline(0)
    %calculate standard deviation
    thisSTD = std(reshape(visobj_history(:,3:end),[],1));
    title(join(["std = " num2str(round(thisSTD))],""))
end

% save
saveas(gcf,'overlap_excite_AOTUmodel.png');

%% vary target noise
% set variables
noiseLevel = 5:5:40;
AOTU019synapse = "inhibitory";
AOTU019overlap = "normal";
nCond = length(noiseLevel);

% initialize
figure
set(gcf,'Position',[100 100 1200 800]);
tiledlayout(nCond/2,6,'TileSpacing','compact')
all_runs = [];
inhibitSTD = [];

% for each overlap condition
for nc = 1:nCond
    % select this condition
    thisNoise = noiseLevel(nc);

    % run this simulation condition
    disp(join(["Simulating " num2str(thisNoise) " noise condition..."],""))
    [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(thisNoise,AOTU019synapse,AOTU019overlap,nTime);
    disp("Simulating complete.")
    % store this simulation condition
    all_runs(:,:,nc) = visobj_history;
    %calculate standard deviation
    inhibitSTD(nc) = std(reshape(visobj_history(:,3:end),[],1));

    % generate condition plot (every other)
    if ~mod(nc,2)
        % plot visual tuning
        nexttile
        plot(visobj_position,AOTU019,"Color","#0072BD")
        hold on
        plot(visobj_position,AOTU025,"Color","#7E2F8E")
        ylabel("firing rate")
        box off,xlim([-180 180]),xticks([-180 0 180])
        % plot example simulation runs
        nexttile([1 3])
        plot(timebase(1,2:end),visobj_history(1:nPlot,2:end),'--o');
        ylim([-180 180]), yticks([-180 0 180]), xticks([])
        yline(0), box off
        title(join([num2str(thisNoise) " noise, " AOTU019synapse],""))
        xlabel('time')
        ylabel('target pos (degrees)')
        % plot polar histogram
        nexttile
        polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a)
        nexttile
        %plot normal histogram
        histogram(visobj_history(:,3:end),'BinWidth',histBinSize,'FaceColor','#77AC30','FaceAlpha',a)
        xlabel('target pos (deg)')
        ylabel('counts')
        xline(0)
        title(join(["std = " num2str(round(inhibitSTD(nc)))],""))
    end
end

% save
saveas(gcf,'noise_AOTUmodel.png');


% repeat with excitatory, wide synapse
AOTU019synapse = "excitatory";
AOTU019overlap = "medium";
nCond = length(noiseLevel);

% initialize
figure
set(gcf,'Position',[100 100 1200 800]);
tiledlayout(nCond/2,6,'TileSpacing','compact')
all_runs = [];
exciteSTD = [];

% for each overlap condition
for nc = 1:nCond
    % select this condition
    thisNoise = noiseLevel(nc);

    % run this simulation condition
    disp(join(["Simulating " num2str(thisNoise) " noise condition w/ excitatory..."],""))
    [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(thisNoise,AOTU019synapse,AOTU019overlap,nTime);
    disp("Simulating complete.")
    % store this simulation condition
    all_runs(:,:,nc) = visobj_history;
    %calculate standard deviation
    exciteSTD(nc) = std(reshape(visobj_history(:,3:end),[],1));

    % generate condition plot (every other)
    if ~mod(nc,2)
        % plot visual tuning
        nexttile
        plot(visobj_position,AOTU019,"Color","#0072BD")
        hold on
        plot(visobj_position,AOTU025,"Color","#7E2F8E")
        ylabel("firing rate")
        box off,xlim([-180 180]),xticks([-180 0 180])
        % plot example simulation runs
        nexttile([1 3])
        plot(timebase(1,2:end),visobj_history(1:nPlot,2:end),'--o');
        ylim([-180 180]), yticks([-180 0 180]), xticks([])
        yline(0), box off
        title(join([num2str(thisNoise) " noise, " AOTU019synapse],""))
        xlabel('time')
        ylabel('target pos (degrees)')
        % plot polar histogram
        nexttile
        polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a)
        nexttile
        %plot normal histogram
        histogram(visobj_history(:,3:end),'BinWidth',histBinSize,'FaceColor','#77AC30','FaceAlpha',a)
        xlabel('target pos (deg)')
        ylabel('counts')
        xline(0)
        title(join(["std = " num2str(round(exciteSTD(nc)))],""))
    end
end

% save
saveas(gcf,'noise_excite_AOTUmodel.png');


% plot inhibitory vs excitatory model comparison
% initialize
figure
set(gcf,'Position',[100 100 400 400]);
plot(noiseLevel,inhibitSTD,"-o","Color","#A2142F")
hold on
plot(noiseLevel,exciteSTD,"-o","Color","#77AC30")
xlabel("noise level")
ylabel("simulation std")

% save
saveas(gcf,'noise_compare_AOTUmodel.png');

