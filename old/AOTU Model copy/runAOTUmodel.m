% runAOTUmodel
% main script for running many itterations of a simple AOTU
% model while modifying basic parameters such as synapse weight, synapse
% sign, noise, open-loop v closed-loop, etc.
%
% 06/20/2024 - MC created
%

%% initialize variables
clear

% set folder path info
thisFile = matlab.desktop.editor.getActiveFilename;
[filePath,~,~] = fileparts(thisFile);
trialsFolder = [filePath '/trials'];
vectorFolder = [trialsFolder '/vectors'];
summaryFolder = [filePath '/summary'];
% fetch receptive field data
cd(filePath)
load("Pursuit_RFs.mat");
nIn = size(Pursuit_RFs,2);

% process/simplify receptive field data
% fetch all cell types
celltypes = Pursuit_RFs.Properties.VariableNames;
% fetch all RFs that are not AOTU019 and AOTU025
idxOther = ~ismember(celltypes,{'AOTU019','AOTU025'});
rfOther = Pursuit_RFs(:,idxOther);
% CB2070 projects contralaterally, flip RF before summing R/L visual inputs
idxCB2070 = ismember(rfOther.Properties.VariableNames,{'CB20701','CB20702'});
rfOther(:,idxCB2070) = flip(rfOther(:,idxCB2070),1);
% sum
rfOther = sum(rfOther,2);
% (optional) lightly smooth combined RF
%gwin = 15;
%rfOther = smoothdata(rfOther,"gaussian",gwin);
% generate RF lookup table for AOTU019, AOTU025, and the sum of all other NOIs
rfAll = [rfOther Pursuit_RFs(:,ismember(celltypes,{'AOTU019','AOTU025'}))];


% set plot variables
% plot - RF variables
rf_ticks = -180:60:180;

% plot - run variables
nEx = 5; % number of example traces

% plot - velocity parameters
vMax = 75;
vBinSize = 5;
vEdges = -vMax+vBinSize/2:vBinSize:vMax;
vCentr = vEdges(1:end-1)+vBinSize/2;

% plot - heading/position parameters
nPolarBins = 20; %polar histogram number of bins
hdBinSize = 10; %normal histogram size of bins

% general plot settings
plotDim = [100 100 1800 800]; %dimensions
celltypes_main = rfAll.Properties.VariableNames;
rfCol = {"#77AC30";"#0072BD";"#7E2F8E"}; %colors for RF curves (sum, 19, 25)
vHistRange = [0 0.2]; %histogram range
hHistRange = [0 0.4]; %histogram range
a = 0.4; %opacity

close all


%% calibrate DNa non-linearity

% initialize
close all
model_settings
thisNoise = 4;
thisStart = 130;
thisSyn = 'excitatory';
simTime = 25;

% run standard model
[timebase,visobj_history,input_history,~] = aotu_steering_model(rfAll,thisNoise,thisStart,thisSyn,simTime);
% determine current offset
cur_offset = DNa02input(interp1(DNa02output,1:length(DNa02output),0.5,'nearest'));
% determine if DN input/output needs to be offset
max_input = max(input_history,[],'all');
min_input = min(input_history,[],'all');
rec_offset = min_input + (max_input - min_input)/2;

% plot
figure; set(gcf,'Position',[100 100 1800 170]);tiledlayout(1,7,'TileSpacing','compact')
% plot visual tuning
nexttile; hold on
for i = 1:3
    plot(visobj_position,rfAll{:,i},"Color",rfCol{i}) %right
    plot(visobj_position,flip(rfAll{:,i}),"Color",rfCol{i}) %left
end
box off;xlim([-180 180]);xticks(rf_ticks);ylabel('FR');xlabel('target pos')
% plot example simulation runs
nexttile([1 3])
plot(timebase(1,2:end),visobj_history(1:nEx,2:end),':.');
axis tight; ylim([-180 180]), yticks([-180 0 180]), xticks([]); yline(0), box off
xlabel('time'); ylabel('target pos (deg)')
% plot polar histogram
nexttile
polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
% plot normal histogram
nexttile
histogram(visobj_history(:,3:end),'BinWidth',hdBinSize,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
xlim([-180 180]); xticks(rf_ticks); ylim(hHistRange); xline(0); xlabel('target pos (deg)')
% plot input histogram
nexttile;
histogram(input_history,40,'FaceAlpha',a,'Normalization','probability')
xlabel('DN Input'); ylabel('norm(probability)')
xline(cur_offset,'g'); xline(rec_offset,'r')

disp(['Current offset: ' num2str(cur_offset)])
disp(['Recommended offset: ' num2str(rec_offset)])

% save image
sgtitle('Test Run')
cd(trialsFolder)
plotname = 'testrun';
saveas(gcf,join([plotname '.png'],''));
% save vectorized plot
cd(vectorFolder)
set(gcf,'renderer','Painters')
saveas(gcf, join([plotname,'.svg'],''))


%% vary k gain and AOTU019 synapse weight
% initialize variables
% vary steering gain (k)
kLevels = 0:10:200;
nK = length(kLevels);
% vary aotu019 synapse strength
strengthLevel = 0:0.2:1;
nStrength = length(strengthLevel);
% generate colormap
cMap = colormap(cool(nStrength));
cMap2 = colormap(turbo(5));

% select loop parameters
close all
thisRun = 'full';
thisStart = 120;
thisNoise = 5;
simTime = 50;
AOTU019synapse = {"inhibitory",'excitatory','slowinhibitory','slowexcitatory'};

% initialize data storage arrays
runSTD = []; runVar = []; runProb = [];

% (optional) compare excitatory and inhibitory performance
nSyn = length(AOTU019synapse);
for y = 1:nSyn
    thisSyn = AOTU019synapse{y};
    as = thisSyn(1); %for storing plots
    disp(join(["Simulating " thisSyn " synapse condition..."],""))

    % for each synapse strength
    for s = 1:nStrength
        thisStrength = strengthLevel(s);
        thisTuning = rfAll;
        thisTuning.AOTU019 = thisTuning.AOTU019.*thisStrength;
        % new figure
        figure; set(gcf,'Position',plotDim);tiledlayout(6,6,'TileSpacing','compact')

        % for each gain
        for n = 1:nK
            % select this condition
            thisK = kLevels(n);

            % run this simulation condition
            disp(join(["Simulating " num2str(thisK) " steering gain condition..."],""))
            [timebase,visobj_history,rotvel_history] = aotu_steering_model_k(thisTuning,thisNoise,thisK,thisStart,thisSyn,simTime);
            disp("Simulation complete.")
            % estimate time to setpoint (0)
            runCross(s,n,y) = median(timebase(findFirstZeroCrossingMultiple(visobj_history)));
            %calculate circular standard deviation
            runSTD(s,n,y) = median(circ_std(deg2rad(visobj_history(:,:)),[],[],2));
            % calculate circular variance
            runVar(s,n,y) = median(1-circ_var(deg2rad(visobj_history(:,:)),[],[],2));
            % pull probability of being within +/- 5 degrees of setpoint
            [h,e] = histcounts(visobj_history(:,:),'BinWidth',hdBinSize,'Normalization','probability');
            e = e(1:end-1)+hdBinSize/2; %center bins
            runProb(s,n,y) = sum(h(abs(e)<=5));

            % generate plot
            if (mod(n,4))==1
                % plot visual tuning
                nexttile; hold on
                for i = 1:3
                    plot(visobj_position,thisTuning{:,i},"Color",rfCol{i}) %right
                    plot(visobj_position,flip(thisTuning{:,i}),"Color",rfCol{i}) %left
                end
                box off;xlim([-180 180]);xticks(rf_ticks);ylabel('FR');xlabel('target pos')

                % plot example simulation runs
                nexttile([1 3])
                plot(timebase(1,2:end),visobj_history(1:nEx,2:end),':.');
                axis tight; ylim([-180 180]), yticks([-180 0 180]), xticks([]); yline(0), box off
                title(join([num2str(thisK) " gain, " thisSyn],""))
                xlabel('time'); ylabel('target pos (deg)')
                % plot polar histogram
                nexttile
                polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                %plot normal histogram
                nexttile
                histogram(visobj_history(:,3:end),'BinWidth',hdBinSize,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                xlim([-180 180]); xticks(rf_ticks); ylim(hHistRange); xline(0); xlabel('target pos (deg)'); ylabel('norm(probability)')
                title(join(["cstd = " num2str(runSTD(s,n,y))],""))
            end
        end
        % save image
        sgtitle(join([thisRun ' Model: AOTU019 vs Steering Gain'],""))
        cd(trialsFolder)
        plotname = join([thisRun 'AOTU019' thisSyn sprintf('%02d',thisStrength*10) 'v_k'],'_');
        saveas(gcf,join([plotname '.png'],''));
        % save vectorized plot
        cd(vectorFolder)
        set(gcf,'renderer','Painters')
        saveas(gcf, join([plotname,'.svg'],''))
    end
end

% plot summary compmaring model performance
figure; set(gcf,'Position',[100,100,1200,900]); tiledlayout(4,4,'TileSpacing','compact')
for y = 1:nSyn
    nexttile; hold on;
    for s = 1:nStrength
        plot(kLevels,runCross(s,:,y),"-o","Color",cMap(s,:))
    end
    ylim([0 6])
    xlabel("k"); ylabel("Time to Setpoint")
    title(AOTU019synapse{y})
end

for y = 1:nSyn
    nexttile; hold on;
    for s = 1:nStrength
        plot(kLevels,runProb(s,:,y),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("k"); ylabel("Probability(near setpoint)")
end

for y = 1:nSyn
    nexttile; hold on;
    for s = 1:nStrength
        plot(kLevels,runVar(s,:,y),"-o","Color",cMap(s,:))
    end
    ylim([0 1])
    xlabel("k"); ylabel("Median 1-Circular Var")
end

for y = 1:nSyn
    nexttile; hold on;
    for s = 1:nStrength
        plot(kLevels,runSTD(s,:,y),"-o","Color",cMap(s,:))
    end
    ylim([0 1.5])
    xlabel("k"); ylabel("Median Circular STD")
end

% save image
sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
cd(summaryFolder)
plotname = join(['summary_' thisRun '_AOTU019' '_v_k'],'');
saveas(gcf,[plotname '.png']);
% save vectorized plot
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])

% plot summary compmaring model performance together
figure; set(gcf,'Position',[100,100,350,900]); tiledlayout(4,1,'TileSpacing','compact')
nexttile; hold on;
for y = 1:nSyn
    plot(kLevels,runCross(nStrength,:,y),"-o","Color",cMap2(y,:))
    ylim([0 6])
    xlabel("k"); ylabel("Time to Setpoint")
end
legend(AOTU019synapse)
nexttile; hold on;
for y = 1:nSyn
    plot(kLevels,runProb(nStrength,:,y),"-o","Color",cMap2(y,:))
    ylim([0 1])
    xlabel("k"); ylabel("Probability(near setpoint)")
end
nexttile; hold on;
for y = 1:nSyn
    plot(kLevels,runVar(nStrength,:,y),"-o","Color",cMap2(y,:))
    ylim([0 1])
    xlabel("k"); ylabel("Median 1-Circular Var")
end
nexttile; hold on;
for y = 1:nSyn
    plot(kLevels,runSTD(nStrength,:,y),"-o","Color",cMap2(y,:))
    ylim([0 1.5])
    xlabel("k"); ylabel("Median Circular STD")
end

% save image
sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
cd(summaryFolder)
plotname = join(['summary_' thisRun '_AOTU019' '_v_kc'],'');
saveas(gcf,[plotname '.png']);
% save vectorized plot
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])

disp('ALL RUNS COMPLETE.')


%% run oscillatory target model varying AOTU019 strength
close all

% set variables
simTime = 975; %time
noiseLevel = 2;
AOTU019synapse = {"inhibitory",'excitatory','slowinhibitory','slowexcitatory'};
nSyn = length(AOTU019synapse);
AOTU019strength = 0:0.2:1;
nStrength = length(AOTU019strength);

% initialize
velDistribution = [];
cMap = colormap(cool(nStrength)); %colormap
cMap2 = colormap(turbo(5));

for y = 1:nSyn
    % select this condition
    thisSynapse = AOTU019synapse{y};
    disp(['BATCH START: AOTU019 = ' thisSynapse])

    % initialize plot
    figure; set(gcf,'Position',plotDim);
    tiledlayout(nStrength,6,'TileSpacing','compact')

    % for each overlap condition
    for nc = 1:nStrength
        % select this condition
        thisStrength = AOTU019strength(nc);
        thisTuning = rfAll;
        thisTuning.AOTU019 = thisTuning.AOTU019.*thisStrength;

        % run this simulation condition
        disp(join(["Simulating " thisStrength "x weight synapse condition..."],""))
        [timebase,visobj_history,rotvel_history] = aotu_steering_oscmodel(thisTuning,noiseLevel,thisSynapse,simTime);
        disp("Simulation complete.")

        % generate condition plot
        % plot visual tuning
        nexttile; hold on
        for i = 1:3
            plot(visobj_position,thisTuning{:,i},"Color",rfCol{i}) %right
            plot(visobj_position,flip(thisTuning{:,i}),"Color",rfCol{i}) %left
        end
        box off;xlim([-180 180]);xticks(rf_ticks);ylabel('FR');xlabel('target pos')

        % plot example simulation runs
        nexttile([1 4])
        plot(timebase,rotvel_history(1:nEx,:),'-');
        axis tight; ylim([-vMax vMax]), xticks([])
        yline(0), box off
        title(join([thisStrength "x weight"],""))
        xlabel('time'); ylabel('velocity')

        %plot normal histogram of velocity
        nexttile
        h = histogram(rotvel_history,'BinEdges',vEdges,'FaceColor','#0072BD','FaceAlpha',a,'Normalization','probability');
        ylim(vHistRange);xlim([-vMax vMax]); xline(0); xlabel('velocity'); ylabel('prob')

        % store
        velDistribution(nc,:,y) = h.Values;
    end
    sgtitle(join(['Oscillatory Model: AOTU019 ' thisSynapse ' v Strength']))
    % save image
    cd(trialsFolder)
    plotname = join(['oscstim_AOTU19',thisSynapse(1)],'_');
    saveas(gcf,join([plotname '.png'],''));
    % save vectorized plot
    cd(vectorFolder)
    set(gcf,'renderer','Painters')
    saveas(gcf, join([plotname '.svg'],''))
end

% plot velocity distribution summary
figure; set(gcf,'Position',[100 100 600 600]);
tiledlayout(2,2,'TileSpacing','compact')
for y = 1:nSyn
    nexttile; hold on
    for s = 1:nStrength
        plot(vCentr,velDistribution(s,:,y),"Color",cMap(s,:),'LineWidth',1.5)
    end
    xlim([-vMax vMax]); xline(0); xlabel('angular velocity'); ylabel('probability')
    title(AOTU019synapse{y})
end
% save image
cd(summaryFolder)
plotname = 'summary_AOTU019_oscstim';
saveas(gcf,[plotname '.png']);
% save vectorized plot
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])

% plot velocity distribution summary together
figure; set(gcf,'Position',[100 100 900 400]);
tiledlayout(1,2,'TileSpacing','compact')
nexttile; hold on
for  y = 1:nSyn
    plot(vCentr,velDistribution(nStrength,:,y),"Color",cMap2(y,:),'LineWidth',1.5)
end
xlim([-vMax vMax]); xline(0); xlabel('angular velocity'); ylabel('probability')
title('1X'); legend(AOTU019synapse)
nexttile; hold on
for  y = 1:nSyn
    plot(vCentr,velDistribution(2,:,y),"Color",cMap2(y,:),'LineWidth',1.5)
end
xlim([-vMax vMax]); xline(0); xlabel('angular velocity'); ylabel('probability')
title('0.2X');
% save image
cd(summaryFolder)
plotname = 'summary_AOTU019_oscstimc';
saveas(gcf,[plotname '.png']);
% save vectorized plot
set(gcf,'renderer','Painters')
saveas(gcf, [plotname '.svg'])


%% vary target noise and AOTU019 synapse weight
% initialize variables
% vary input noise level
noiseLevel = 0:2:40;
nNoise = length(noiseLevel);
% vary aotu019 synapse strength
strengthLevel = 0:0.2:1;
nStrength = length(strengthLevel);
% generate colormap
cMap = colormap(cool(nStrength));

% set loop parameters
close all
runOpt = {'correction';'stabilize'};
startOpt = [120, 0];
timeOpt = [25, 25];
AOTU019synapse = {"inhibitory",'excitatory','slowinhibitory','slowexcitatory'};
nSyn = length(AOTU019synapse);

% using parallel processing...
% run for full, correction only, and setpoint only
for r = 1:2
    % select loop parameters
    thisRun = runOpt{r};
    thisStart = startOpt(r);
    simTime = timeOpt(r);

    % initialize data storage arrays
    runSTD = []; runVar = []; runProb = [];

    % (optional) compare excitatory and inhibitory performance
    for y = 1:nSyn
        thisSyn = AOTU019synapse{y};
        disp(join(["Simulating " thisSyn " synapse condition..."],""))

        % for each synapse strength
        for s = 1:nStrength
            thisStrength = strengthLevel(s);
            thisTuning = rfAll;
            thisTuning.AOTU019 = thisTuning.AOTU019.*thisStrength;
            % new figure
            figure; set(gcf,'Position',plotDim);tiledlayout(5,6,'TileSpacing','compact')

            % for each noise level
            for n = 1:nNoise
                % select this condition
                thisNoise = noiseLevel(n);

                % run this simulation condition
                disp(join(["Simulating " num2str(thisNoise) " noise condition..."],""))
                [timebase,visobj_history,rotvel_history] = aotu_steering_model(thisTuning,thisNoise,thisStart,thisSyn,simTime);
                disp("Simulation complete.")
                % estimate time to setpoint (0)
                runCross(s,n,y) = median(timebase(findFirstZeroCrossingMultiple(visobj_history)));
                %calculate circular standard deviation
                runSTD(s,n,y) = median(circ_std(deg2rad(visobj_history(:,3:end)),[],[],2));
                % calculate circular variance
                runVar(s,n,y) = median(1-circ_var(deg2rad(visobj_history(:,3:end)),[],[],2));
                % pull probability of being within +/- 5 degrees of setpoint
                [h,e] = histcounts(visobj_history(:,3:end),'BinWidth',hdBinSize,'Normalization','probability');
                e = e(1:end-1)+hdBinSize/2; %center bins
                runProb(s,n,y) = sum(h(abs(e)<=5));

                % generate plot
                if (mod(n,5))==1
                    % plot visual tuning
                    nexttile; hold on
                    for i = 1:3
                        plot(visobj_position,thisTuning{:,i},"Color",rfCol{i}) %right
                        plot(visobj_position,flip(thisTuning{:,i}),"Color",rfCol{i}) %left
                    end
                    box off;xlim([-180 180]);xticks(rf_ticks);ylabel('FR');xlabel('target pos')

                    % plot example simulation runs
                    nexttile([1 3])
                    plot(timebase(1,2:end),visobj_history(1:nEx,2:end),':.');
                    axis tight; ylim([-180 180]), yticks([-180 0 180]), xticks([]); yline(0), box off
                    title(join([num2str(thisNoise) " noise, " thisSyn],""))
                    xlabel('time'); ylabel('target pos (deg)')
                    % plot polar histogram
                    nexttile
                    polarhistogram(deg2rad(visobj_history(:,3:end)),nPolarBins,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                    %plot normal histogram
                    nexttile
                    histogram(visobj_history(:,3:end),'BinWidth',hdBinSize,'FaceColor','#77AC30','FaceAlpha',a,'Normalization','probability')
                    xlim([-180 180]); xticks(rf_ticks); ylim(hHistRange); xline(0); xlabel('target pos (deg)'); ylabel('norm(probability)')
                    title(join(["cstd = " num2str(runSTD(s,n,y))],""))
                end
            end
            % save image
            sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
            cd(trialsFolder)
            plotname = join([thisRun 'AOTU019' thisSyn sprintf('%02d',thisStrength*10) 'v_noise'],'_');
            saveas(gcf,join([plotname '.png'],''));
            % save vectorized plot
            cd(vectorFolder)
            set(gcf,'renderer','Painters')
            saveas(gcf, join([plotname,'.svg'],''))
        end
    end

    % plot summary compmaring model performance
    if r==1
        nt = 4;
        thisPlotdim = [100,100,900,900];
    else
        nt = 3;
        thisPlotdim = [100,100,900,900];
    end
    figure; set(gcf,'Position',thisPlotdim); tiledlayout(nt,nSyn,'TileSpacing','compact')
    if r==1
        for y = 1:nSyn
            nexttile; hold on;
            for s = 1:nStrength
                plot(noiseLevel,runCross(s,:,y),"-o","Color",cMap(s,:))
            end
            ylim([0 1.5])
            xlabel("k"); ylabel("Time to Setpoint")
            title(AOTU019synapse{y})
        end
    end
    for y = 1:nSyn
        nexttile; hold on;
        for s = 1:nStrength
            plot(noiseLevel,runProb(s,:,y),"-o","Color",cMap(s,:))
        end
        ylim([0 1])
        xlabel("noise level"); ylabel("Probability(near setpoint)")
        title(AOTU019synapse{y})
    end

    for y = 1:nSyn
        nexttile; hold on;
        for s = 1:nStrength
            plot(noiseLevel,runVar(s,:,y),"-o","Color",cMap(s,:))
        end
        ylim([0 1])
        xlabel("noise level"); ylabel("Median 1-Circular Var")
    end

    for y = 1:nSyn
        nexttile; hold on;
        for s = 1:nStrength
            plot(noiseLevel,runSTD(s,:,y),"-o","Color",cMap(s,:))
        end
        ylim([0 1.5])
        xlabel("noise level"); ylabel("Median Circular STD")
    end

    % save image
    sgtitle(join([thisRun ' Model: AOTU019 vs Stimulus Noise'],""))
    cd(summaryFolder)
    plotname = join(['summary_' thisRun '_AOTU019' '_v_noise'],'');
    saveas(gcf,[plotname '.png']);
    % save vectorized plot
    set(gcf,'renderer','Painters')
    saveas(gcf, [plotname '.svg'])

end
disp('ALL RUNS COMPLETE.')

