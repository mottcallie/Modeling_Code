% runAOTUmodel
% main script for running many itterations of a simple AOTU019 + AOTU025
% model while modifying basic parameters such as synapse weight, synapse
% sign, noise, open-loop v closed-loop, etc.
%
% 05/01/2024 - MC/RW created
% 05/06/2024 - MC added circular statistics
% 05/20/2024 - MC added oscillating target model
%

%% initialize variables
clear
close all

% set folder path info
thisFile = matlab.desktop.editor.getActiveFilename;
[filePath,~,~] = fileparts(thisFile);
plotFolder = [filePath '/plots'];
vectorFolder = [plotFolder '/vectors'];
plotFolder2 = [filePath '/plots_l'];
oscFolder = [filePath '/oscillating'];
vectorFolder2 = [plotFolder2 '/vectors'];
summaryFolder = [filePath '/summary'];
cd([filePath '/plots'])

% plot - RF variables
visobj_position = linspace(-180, 180, 361); %obj pos
% plot - run variables
nEx = 5; % number of example traces
% plot - histogram variables
a = 0.4; %opacity
nPolarBins = 20; %polar histogram number of bins
hdBinSize = 10; %normal histogram size of bins
velBinSize = 5; %normal histogram size of bins
plotDim = [100 100 1800 800]; %dimensions
velBins = 1;
velEdges = -100+velBins/2:velBins:100;
velCentr = velEdges(1:end-1)+velBins/2;


%% vary AOTU019 synapse speed
% set variables
noiseLevel = 10;
thisSyn = {"inhibitory","slow"};
AOTU019strength = 1;
AOTU019overlap = "normal";
DNa02trans = "non";
nCond = length(thisSyn);
thisTime = 80; %time of simulation run
thisStart = 120; %+/-deg

% initialize
figure
set(gcf,'Position',plotDim);
tiledlayout(nCond,6,'TileSpacing','compact')
all_runs = [];

% for each overlap condition
for nc = 1:nCond
    % select this condition
    thisSynapse = thisSyn{nc};

    % run this simulation condition
    disp(join(["Simulating " thisSynapse " synapse condition..."],""))
    [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(noiseLevel,thisStart,thisSynapse,AOTU019strength,AOTU019overlap,DNa02trans,thisTime);
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
    plot(timebase(1,2:end),visobj_history(1:nEx,2:end),':o');
    ylim([-180 180]), yticks([-180 0 180]), xticks([])
    yline(0), box off
    title(join([thisSynapse " synapse"],""))
    xlabel('time')
    ylabel('target pos (degrees)')
    % plot polar histogram
    nexttile
    polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
    %plot normal histogram
    nexttile
    histogram(visobj_history(:,3:end),'BinWidth',hdBinSize,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
    xlabel('target pos (deg)')
    ylabel('norm(probability)')
    xline(0)
    %calculate standard deviation
    thisSTD = median(circ_std(deg2rad(visobj_history(:,3:end)),[],[],2));
    title(join(["cstd = " num2str(thisSTD)],""))
end
sgtitle('Model: AOTU019 Properties')

% save image
cd(plotFolder)
plotname = 'full_AOTU_synmodel';
saveas(gcf,[plotname '.png']);
% save vectorized plot
cd(vectorFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])



%% vary target noise and AOTU019 synapse weight
close all
% initialize variables
% vary input noise level
noiseLevel = 0:2:50;
nNoise = length(noiseLevel);
% vary aotu019 synapse strength
strengthLevel = 0:0.2:1;
nStrength = length(strengthLevel);
AOTU019overlap = "normal";
DNa02trans = "non";
% generate colormap
cMap = colormap(cool(nStrength));

% set parallel loop parameters
runOpt = {'full';'correct';'setpoint'};
startOpt = [120, 120, 0];
timeOpt = [80, 25, 25];

synapseOpt = {'inhibitory';'excitatory'};

% using parallel processing...
% run for full, correction only, and setpoint only
for r = 1:3
    % select loop parameters
    thisRun = runOpt{r};
    thisStart = startOpt(r);
    thisTime = timeOpt(r);

    % initialize data storage arrays
    runSTD = []; runVar = []; runProb = [];

    % compare excitatory and inhibitory performance
    for y = 1:2
        thisSyn = synapseOpt{y};
        as = thisSyn(1); %for storing plots
        disp(join(["Simulating " thisSyn " synapse condition..."],""))

        % for each synapse strength
        for s = 1:nStrength
            thisStrength = strengthLevel(s);
            % new figure
            figure; set(gcf,'Position',plotDim);tiledlayout(5,6,'TileSpacing','compact')

            % for each noise level
            for n = 1:nNoise
                % select this condition
                thisNoise = noiseLevel(n);

                % run this simulation condition
                disp(join(["Simulating " num2str(thisNoise) " noise condition..."],""))
                [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(thisNoise,thisStart,thisSyn,thisStrength,AOTU019overlap,DNa02trans,thisTime);
                disp("Simulation complete.")
                %calculate circular standard deviation
                runSTD(s,n,y) = median(circ_std(deg2rad(visobj_history(:,3:end)),[],[],2));
                % calculate circular variance
                runVar(s,n,y) = median(1-circ_var(deg2rad(visobj_history(:,3:end)),[],[],2));
                % pull probability of being within +/- 5 degrees of setpoint
                [h,e] = histcounts(visobj_history(:,3:end),'BinWidth',hdBinSize,'Normalization','probability');
                e = e(1:end-1)+hdBinSize/2; %center bins
                runProb(s,n,y) = sum(h(abs(e)<=5));

                % generate plot
                if (mod(n,6))==1
                    % plot visual tuning
                    nexttile
                    hold on; plot(visobj_position,AOTU019,"Color","#0072BD")
                    plot(visobj_position,AOTU025,"Color","#7E2F8E")
                    ylabel("firing rate"); box off,xlim([-180 180]),xticks([-180 0 180])
                    % plot example simulation runs
                    nexttile([1 3])
                    plot(timebase(1,2:end),visobj_history(1:nEx,2:end),':.');
                    ylim([-180 180]), yticks([-180 0 180]), xticks([]); yline(0), box off
                    title(join([num2str(thisNoise) " noise, " thisSyn],""))
                    xlabel('time'); ylabel('target pos (degrees)')
                    % plot polar histogram
                    nexttile
                    polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                    %plot normal histogram
                    nexttile
                    histogram(visobj_history(:,3:end),'BinWidth',hdBinSize,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                    xlim([-180 180]); xline(0); xlabel('target pos (deg)'); ylabel('norm(probability)')
                    title(join(["cstd = " num2str(runSTD(s,n,y))],""))
                end
            end
            % save image
            sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
            cd(plotFolder)
            plotname = join([thisRun '_AOTU' as '_' sprintf('%02d',thisStrength*10) '_noisemodel'],'');
            saveas(gcf,[plotname '.png']);
            % save vectorized plot
            cd(vectorFolder)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
        end
    end
    % store summary
    sumname = join([thisRun '_summary_data'],'');
    parsave([sumname '.mat'],runProb,runVar,runSTD);

    % plot summary compmaring model performance
    figure; set(gcf,'Position',[100,100,800,900]); tiledlayout(3,2,'TileSpacing','compact')

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runProb(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Probability(near setpoint)")
    title('AOTU019 Inhibitory')
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runProb(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Probability(near setpoint)")
    title('AOTU019 Excitatory')

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runVar(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Median 1-Circular Var")
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runVar(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Median 1-Circular Var")

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runSTD(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1.5])
    xlabel("noise level"); ylabel("Median Circular STD")
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runSTD(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1.5])
    xlabel("noise level"); ylabel("Median Circular STD")

    % save image
    sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
    cd(summaryFolder)
    plotname = join([thisRun '_summary_AOTU' '_noisemodel'],'');
    saveas(gcf,[plotname '.png']);
    % save vectorized plot
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
end
disp('ALL RUNS COMPLETE.')

%% vary target noise and synapse weight with DNa02 linear
% initialize variables
% vary input noise level
noiseLevel = 0:2:50;
nNoise = length(noiseLevel);
% vary aotu019 synapse strength
strengthLevel = 0:0.2:1;
nStrength = length(strengthLevel);
AOTU019overlap = "normal";
DNa02trans = "linear";
% generate colormap
cMap = colormap(cool(nStrength));

% set parallel loop parameters
runOpt = {'full';'correct';'setpoint'};
startOpt = [120, 120, 0];
timeOpt = [80, 25, 25];

synapseOpt = {'inhibitory';'excitatory'};

% using parallel processing...
% run for full, correction only, and setpoint only
parfor r = 1:3
    % select loop parameters
    thisRun = runOpt{r};
    thisStart = startOpt(r);
    thisTime = timeOpt(r);

    % initialize data storage arrays
    runSTD = []; runVar = []; runProb = [];

    % compare excitatory and inhibitory performance
    for y = 1:2
        thisSyn = synapseOpt{y};
        as = thisSyn(1); %for storing plots
        disp(join(["Simulating " thisSyn " synapse condition..."],""))

        % for each synapse strength
        for s = 1:nStrength
            thisStrength = strengthLevel(s);
            % new figure
            figure; set(gcf,'Position',plotDim);tiledlayout(5,6,'TileSpacing','compact')

            % for each noise level
            for n = 1:nNoise
                % select this condition
                thisNoise = noiseLevel(n);

                % run this simulation condition
                disp(join(["Simulating " num2str(thisNoise) " noise condition..."],""))
                [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(thisNoise,thisStart,thisSyn,thisStrength,AOTU019overlap,DNa02trans,thisTime);
                disp("Simulation complete.")
                %calculate circular standard deviation
                runSTD(s,n,y) = median(circ_std(deg2rad(visobj_history(:,3:end)),[],[],2));
                % calculate circular variance
                runVar(s,n,y) = median(1-circ_var(deg2rad(visobj_history(:,3:end)),[],[],2));
                % pull probability of being within +/- 5 degrees of setpoint
                [h,e] = histcounts(visobj_history(:,3:end),'BinWidth',hdBinSize,'Normalization','probability');
                e = e(1:end-1)+hdBinSize/2; %center bins
                runProb(s,n,y) = sum(h(abs(e)<=5));

                % generate plot
                if (mod(n,6))==1
                    % plot visual tuning
                    nexttile
                    hold on; plot(visobj_position,AOTU019,"Color","#0072BD")
                    plot(visobj_position,AOTU025,"Color","#7E2F8E")
                    ylabel("firing rate"); box off,xlim([-180 180]),xticks([-180 0 180])
                    % plot example simulation runs
                    nexttile([1 3])
                    plot(timebase(1,2:end),visobj_history(1:nEx,2:end),':.');
                    ylim([-180 180]), yticks([-180 0 180]), xticks([]); yline(0), box off
                    title(join([num2str(thisNoise) " noise, " thisSyn],""))
                    xlabel('time'); ylabel('target pos (degrees)')
                    % plot polar histogram
                    nexttile
                    polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                    %plot normal histogram
                    nexttile
                    histogram(visobj_history(:,3:end),'BinWidth',hdBinSize,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                    xlim([-180 180]); xline(0); xlabel('target pos (deg)'); ylabel('norm(probability)')
                    title(join(["cstd = " num2str(runSTD(s,n,y))],""))
                end
            end
            % save image
            sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise (a02 linear)'],""))
            cd(plotFolder2)
            plotname = join([thisRun '_AOTU' as '_' sprintf('%02d',thisStrength*10) '_lin_noisemodel'],'');
            saveas(gcf,[plotname '.png']);
            % save vectorized plot
            cd(vectorFolder)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
        end
    end
    % store summary
    sumname = join([thisRun '_summary_data'],'');
    parsave([sumname '.mat'],runProb,runVar,runSTD);

    % plot summary compmaring model performance
    figure; set(gcf,'Position',[100,100,800,900]); tiledlayout(3,2,'TileSpacing','compact')

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runProb(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Probability(near setpoint)")
    title('AOTU019 Inhibitory')
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runProb(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Probability(near setpoint)")
    title('AOTU019 Excitatory')

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runVar(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Median 1-Circular Var")
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runVar(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Median 1-Circular Var")

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runSTD(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1.5])
    xlabel("noise level"); ylabel("Median Circular STD")
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runSTD(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1.5])
    xlabel("noise level"); ylabel("Median Circular STD")

    % save image
    sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
    cd(summaryFolder)
    plotname = join([thisRun '_summary_AOTU' '_lin_noisemodel'],'');
    saveas(gcf,[plotname '.png']);
    % save vectorized plot
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
end
disp('ALL RUNS COMPLETE.')

%% vary target noise and AOTU019 synapse weight

% initialize variables
% vary input noise level
noiseLevel = 0:2:50;
nNoise = length(noiseLevel);
% vary aotu019 synapse strength
strengthLevel = 0:0.2:1;
nStrength = length(strengthLevel);
AOTU019overlap = "normal";
DNa02trans = "non";
% generate colormap
cMap = colormap(cool(nStrength));

% set parallel loop parameters
runOpt = {'full';'correct';'setpoint'};
startOpt = [120, 120, 0];
timeOpt = [80, 25, 25];

synapseOpt = {'inhibitory';'excitatory'};

% using parallel processing...
% run for full, correction only, and setpoint only
parfor r = 1:3
    % select loop parameters
    thisRun = runOpt{r};
    thisStart = startOpt(r);
    thisTime = timeOpt(r);

    % initialize data storage arrays
    runSTD = []; runVar = []; runProb = [];

    % compare excitatory and inhibitory performance
    for y = 1:2
        thisSyn = synapseOpt{y};
        as = thisSyn(1); %for storing plots
        disp(join(["Simulating " thisSyn " synapse condition..."],""))

        % for each synapse strength
        for s = 1:nStrength
            thisStrength = strengthLevel(s);
            % new figure
            figure; set(gcf,'Position',plotDim);tiledlayout(5,6,'TileSpacing','compact')

            % for each noise level
            for n = 1:nNoise
                % select this condition
                thisNoise = noiseLevel(n);

                % run this simulation condition
                disp(join(["Simulating " num2str(thisNoise) " noise condition..."],""))
                [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model(thisNoise,thisStart,thisSyn,thisStrength,AOTU019overlap,DNa02trans,thisTime);
                disp("Simulation complete.")
                %calculate circular standard deviation
                runSTD(s,n,y) = median(circ_std(deg2rad(visobj_history(:,3:end)),[],[],2));
                % calculate circular variance
                runVar(s,n,y) = median(1-circ_var(deg2rad(visobj_history(:,3:end)),[],[],2));
                % pull probability of being within +/- 5 degrees of setpoint
                [h,e] = histcounts(visobj_history(:,3:end),'BinWidth',hdBinSize,'Normalization','probability');
                e = e(1:end-1)+hdBinSize/2; %center bins
                runProb(s,n,y) = sum(h(abs(e)<=5));

                % generate plot
                if (mod(n,6))==1
                    % plot visual tuning
                    nexttile
                    hold on; plot(visobj_position,AOTU019,"Color","#0072BD")
                    plot(visobj_position,AOTU025,"Color","#7E2F8E")
                    ylabel("firing rate"); box off,xlim([-180 180]),xticks([-180 0 180])
                    % plot example simulation runs
                    nexttile([1 3])
                    plot(timebase(1,2:end),visobj_history(1:nEx,2:end),':.');
                    ylim([-180 180]), yticks([-180 0 180]), xticks([]); yline(0), box off
                    title(join([num2str(thisNoise) " noise, " thisSyn],""))
                    xlabel('time'); ylabel('target pos (degrees)')
                    % plot polar histogram
                    nexttile
                    polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                    %plot normal histogram
                    nexttile
                    histogram(visobj_history(:,3:end),'BinWidth',hdBinSize,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                    xlim([-180 180]); xline(0); xlabel('target pos (deg)'); ylabel('norm(probability)')
                    title(join(["cstd = " num2str(runSTD(s,n,y))],""))
                end
            end
            % save image
            sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
            cd(plotFolder)
            plotname = join([thisRun '_AOTU' as '_' sprintf('%02d',thisStrength*10) '_noisemodel'],'');
            saveas(gcf,[plotname '.png']);
            % save vectorized plot
            cd(vectorFolder)
            set(gcf,'renderer','Painters')
            saveas(gcf, [plotname '.svg'])
        end
    end
    % store summary
    sumname = join([thisRun '_summary_data'],'');
    parsave([sumname '.mat'],runProb,runVar,runSTD);

    % plot summary compmaring model performance
    figure; set(gcf,'Position',[100,100,800,900]); tiledlayout(3,2,'TileSpacing','compact')

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runProb(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Probability(near setpoint)")
    title('AOTU019 Inhibitory')
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runProb(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Probability(near setpoint)")
    title('AOTU019 Excitatory')

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runVar(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Median 1-Circular Var")
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runVar(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("noise level"); ylabel("Median 1-Circular Var")

    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runSTD(s,:,1),"-o","Color",cMap(s,:))
    end
    ylim([0 1.5])
    xlabel("noise level"); ylabel("Median Circular STD")
    nexttile; hold on;
    for s = 1:nStrength
        plot(noiseLevel,runSTD(s,:,2),"-o","Color",cMap(s,:))
    end
    ylim([0 1.5])
    xlabel("noise level"); ylabel("Median Circular STD")

    % save image
    sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
    cd(summaryFolder)
    plotname = join([thisRun '_summary_AOTU' '_noisemodel'],'');
    saveas(gcf,[plotname '.png']);
    % save vectorized plot
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])
end
disp('ALL RUNS COMPLETE.')
close all


%% vary target motion frequency

% initialize variables
% vary lowpass
lowPass = 0:1:5;
nLP = length(lowPass);
% vary input noise level
noiseLevel = 0:2:50;
nNoise = length(noiseLevel);
% vary aotu019 synapse strength
strengthLevel = [0 1];
nStrength = length(strengthLevel);
% generate colormap
cMap = colormap(cool(nStrength));

% for each lowpass
close all
figure; set(gcf,'Position',[100 100 1800 400])
runDiff = [];
for r = 1:nLP
    % select loop parameters
    thisLowPass = lowPass(r);
    thisRun = 'setpoint';
    thisStart = 0;
    thisTime = 50;

    % initialize data storage arrays
    runSTD = []; runVar = []; runProb = [];

    % compare inhibitory performance
    thisSyn = 'inhibitory';
    as = thisSyn(1); %for storing plots
    disp(join(["Simulating " thisSyn " synapse condition..."],""))

    % for each synapse strength
    for s = 1:nStrength
        thisStrength = strengthLevel(s);

        % for each noise level
        for n = 1:nNoise
            % select this condition
            thisNoise = noiseLevel(n);

            % run this simulation condition
            disp(join(["Simulating " num2str(thisNoise) " noise condition..."],""))
            [AOTU019,AOTU025,timebase,visobj_history] = aotu_steering_model_lowpass(thisLowPass,thisNoise,thisStart,thisSyn,thisStrength,thisTime);
            disp("Simulation complete.")
            %calculate circular standard deviation
            runSTD(s,n) = median(circ_std(deg2rad(visobj_history(:,3:end)),[],[],2));
            % calculate circular variance
            runVar(s,n) = median(1-circ_var(deg2rad(visobj_history(:,3:end)),[],[],2));
            % pull probability of being within +/- 5 degrees of setpoint
            [h,e] = histcounts(visobj_history(:,3:end),'BinWidth',hdBinSize,'Normalization','probability');
            e = e(1:end-1)+hdBinSize/2; %center bins
            runProb(s,n) = sum(h(abs(e)<=5));
        end

        % plot setpoint probability
        subplot(1,nLP,r); hold on;
        plot(noiseLevel,runProb(s,:),"-o","Color",cMap(s,:))
        axis tight; ylim([0 1])
        xlabel("Noise Level"); ylabel("Probability(near setpoint)");
        title([num2str(thisLowPass) 'Hz Lowpass'])

    end
    % pull difference between strength performances
    runDiff(r,:) = runProb(2,:) - runProb(1,:);
end
% save image
sgtitle('Model: AOTU019 vs Stimulus Frequency')
cd(summaryFolder)
plotname = join(['freq_summary_AOTU' '_noisemodel'],'');
saveas(gcf,[plotname '.png']);
% save vectorized plot
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])

% plot summary
figure; set(gcf,'Position',[100 100 400 400])
% generate colormap
dMap = colormap(flip(summer(nLP)));
for r = 1:nLP
    hold on
    plot(noiseLevel,runDiff(r,:),'-o',"Color",dMap(r,:))
end
axis padded;
legend(num2str(lowPass')); ylabel('Performance Difference'); xlabel ('Noise Level')

% save image
sgtitle('Model: AOTU019 vs Stimulus Frequency')
cd(summaryFolder)
plotname = join(['freq_summary_AOTU' '_noisemodel'],'');
saveas(gcf,[plotname '.png']);
% save vectorized plot
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])

disp('ALL RUNS COMPLETE.')

%% run oscillatory target model varying AOTU019 strength

% set variables
noiseLevel = 2;
AOTU019synapse = "inhibitory";
AOTU019strength = 0:0.2:1;
AOTU025strength = 1;
nStrength = length(AOTU019strength);
thisTime = 400; %time of simulation run
% generate colormap
cMap = colormap(cool(nStrength));

% initialize
figure; set(gcf,'Position',plotDim);
tiledlayout(nStrength,6,'TileSpacing','compact')
velDistribution = [];

% for each overlap condition
for nc = 1:nStrength
    % select this condition
    thisStrength = AOTU019strength(nc);

    % run this simulation condition
    disp(join(["Simulating " thisStrength "x weight synapse condition..."],""))
    [AOTU019,AOTU025,timebase,visobj_history,rotvel_history] = aotu_steering_oscmodel(noiseLevel,AOTU019synapse,thisStrength,AOTU025strength,thisTime);
    disp("Simulating complete.")

    % generate condition plot
    % plot visual tuning
    nexttile
    plot(visobj_position,AOTU019,"Color","#0072BD")
    hold on
    plot(visobj_position,AOTU025,"Color","#7E2F8E")
    ylabel("firing rate")
    box off,xlim([-180 180]),xticks([-180 0 180])
    % plot example simulation runs
    nexttile([1 4])
    plot(timebase,rotvel_history(1:nEx,:),'-');
    axis padded; ylim([-60 60]), yticks([-50 0 50]), xticks([])
    yline(0), box off
    title(join([thisStrength "x weight"],""))
    xlabel('time')
    ylabel('rotational velocity (deg/s)')
    %plot normal histogram of velocity
    nexttile
    h = histogram(rotvel_history,'BinEdges',velEdges,'FaceColor','#0072BD','FaceAlpha',a,'Normalization','probability');
    ylim([0 0.1]);xlim([-60 60]); xline(0); xlabel('rotational velocity (deg/s)'); ylabel('norm(probability)')
    % store
    velDistribution(nc,:) = h.Values;
end
sgtitle('Oscillatory Model: AOTU019 Strength')

% save image
cd(oscFolder)
plotname = 'osc_AOTU_model';
saveas(gcf,[plotname '.png']);
% save vectorized plot
cd(oscFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])

% plot velocity distribution summary
figure; set(gcf,'Position',[100 100 400 400]); hold on
for s = 1:nStrength
    plot(velCentr,velDistribution(s,:),"Color",cMap(s,:),'LineWidth',1.5)
end
xlim([-80 80]); xline(0); xlabel('rotational velocity (deg/s)'); ylabel('norm(probability)')

% save image
cd(oscFolder)
plotname = 'osc_AOTU_modelsum';
saveas(gcf,[plotname '.png']);
% save vectorized plot
cd(oscFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])

%% run oscillatory target model varying AOTU025

% set variables
noiseLevel = 5;
AOTU019synapse = "inhibitory";
AOTU019strength = 1;
AOTU025strength = 0:0.2:1;
nStrength = length(AOTU025strength);
thisTime = 400; %time of simulation run
% generate colormap
cMap = colormap(cool(nStrength));

% initialize
figure; set(gcf,'Position',plotDim);
tiledlayout(nStrength,6,'TileSpacing','compact')
velDistribution = [];

% for each overlap condition
for nc = 1:nStrength
    % select this condition
    thisStrength = AOTU025strength(nc);

    % run this simulation condition
    disp(join(["Simulating " thisStrength "x weight synapse condition..."],""))
    [AOTU019,AOTU025,timebase,visobj_history,rotvel_history] = aotu_steering_oscmodel(noiseLevel,AOTU019synapse,AOTU019strength,thisStrength,thisTime);
    disp("Simulating complete.")

    % generate condition plot
    % plot visual tuning
    nexttile
    plot(visobj_position,AOTU019,"Color","#0072BD")
    hold on
    plot(visobj_position,AOTU025,"Color","#7E2F8E")
    ylabel("firing rate")
    box off,xlim([-180 180]),xticks([-180 0 180])
    % plot example simulation runs
    nexttile([1 4])
    plot(timebase,rotvel_history(1:nEx,:),'-');
    axis padded; ylim([-60 60]), yticks([-50 0 50]), xticks([])
    yline(0), box off
    title(join([thisStrength "x weight"],""))
    xlabel('time')
    ylabel('rotational velocity (deg/s)')
    %plot normal histogram of velocity
    nexttile
    h = histogram(rotvel_history,'BinEdges',velEdges,'FaceColor','#0072BD','FaceAlpha',a,'Normalization','probability');
    ylim([0 0.1]);xlim([-60 60]); xline(0); xlabel('rotational velocity (deg/s)'); ylabel('norm(probability)')
    % store
    velDistribution(nc,:) = h.Values;
end
sgtitle('Oscillatory Model: AOTU025 Strength')

% save image
cd(oscFolder)
plotname = 'osc_AOTU25_model';
saveas(gcf,[plotname '.png']);
% save vectorized plot
cd(oscFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])

% plot velocity distribution summary
figure; set(gcf,'Position',[100 100 400 400]); hold on
for s = 1:nStrength
    plot(velCentr,velDistribution(s,:),"Color",cMap(s,:),'LineWidth',1.5)
end
xlim([-80 80]); xline(0); xlabel('rotational velocity (deg/s)'); ylabel('norm(probability)')

% save image
cd(oscFolder)
plotname = 'osc_AOTU25_modelsum';
saveas(gcf,[plotname '.png']);
% save vectorized plot
cd(oscFolder)
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])