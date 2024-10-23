% model_pursuit
% script for generating a simple linear model for pursuit beahvior
%
% 04/18/24 MC created
% 05/06/24 MC adjusted fit to 80% of data, measure against 20%, repeat
%

%% initialize

% clean up
clear
close all

% fetch processing settings
settings = processSettings();

% set data folders
thisFile = matlab.desktop.editor.getActiveFilename;
thisFolder = extractBefore(thisFile,'\model_pursuit.m');
plotFolder = [thisFolder '\plots'];
pursuitFolder = "E:\AOTU019 Oscillate\interpolated";
darknssFolder = "E:\AOTU019 Background P1\interpolated";

% load file/fly names from pursuit data
cd(pursuitFolder)
pursuitFiles = dir('*int.mat');
obj_xcorr.Obj = 0;

for n = 1:length(pursuitFiles)
    flyNames{n,1} = extractBefore(pursuitFiles(n).name,'_int.mat');
end

% load corresponding background data
cd(darknssFolder)
darknssFiles = dir('*pulse.mat');
load('trackLR.mat')
trackLR = flip(trackLR);
for n = length(darknssFiles):-1:1
    if ~contains(darknssFiles(n).name,flyNames)
        darknssFiles(n) = [];
        trackLR(n) = [];
    end
end

% set variables
restThreshold = 1; %mm/s
pursThreshold = 3; %mm/s
speedSelect = 5;

cutoffD = 25; %s, min time for each behavior parameter to be included

nFliesP = size(pursuitFiles,1);
nFliesD = size(darknssFiles,1);

% plot variables
trial_grey = [0.5 0.5 0.5];
v_names = {'Forward', 'Ipsi Angular','Ipsi Sideways'};
v_colors = {'#D95319';'#0072BD';'#7E2F8E';'#A2142F'};
p_names = {'Ipsi Target Position at Rest';'Ipsi Target Direction at Rest'};
p_colors = {'#77AC30';'#4DBEEE'};

%% generate linear fit
disp('Fitting sequential linear models...')
close all

% initialize data storage arrays
allVelocity = cell(nFliesD,3);
allPanelps = cell(nFliesD,2);
allSpikert = cell(nFliesD,2);
predictor_t = nan(4,nFliesD);
fasR2 = nan(3,nFliesD);
asfR2 = nan(3,nFliesD);
sfaR2 = nan(3,nFliesD);

pdR2 = nan(2,nFliesD);
dpR2 = nan(2,nFliesD);

fullR2_fapd = nan(4,nFliesD);
fullR2_pdfa = nan(4,nFliesD);
nTrial = zeros(nFliesD,2);

for nf = 1:nFliesD
    disp(['Fitting fly ' num2str(nf) '/' num2str(nFliesD)])
    % first, fit directional velocity from darkness data
    thisFly = extractBefore(darknssFiles(nf).name,'_pulse.mat');
    % fetch darkness dataset
    thisFile = fullfile(darknssFolder,darknssFiles(nf).name);
    dfD = load(thisFile);
    
    % pull response data
    spikert = dfD.int_spikert(:,:,1);
    % pull predictor data and lag shift
    t = dfD.int_time;
    forward = lagShift(dfD.int_forward(:,:,1),t,settings.fwdLag);
    angular = lagShift(dfD.int_angular(:,:,1),t,settings.angLag);
    sideway = lagShift(dfD.int_sideway(:,:,1),t,settings.sidLag);

    % omit rests
    forward(forward<0.25) = nan;
    angular(angular<5) = nan;
    sideway(sideway<0.1) = nan;
    % for turning pull IPSI only
    angular(angular<0) = nan;
    sideway(sideway<0) = nan;

    % fit linear models for directional velocity
    predictorOrder1(1,:) = {'+Forward','+Angular','+Sideways'};
    [mdlFAS] = sequential_lm(spikert,forward,angular,sideway);
    predictorOrder1(2,:) = {'+Angular','+Sideways','+Forward'};
    [mdlASF] = sequential_lm(spikert,angular,sideway,forward);
    predictorOrder1(3,:) = {'+Sideways','+Forward','+Angular'};
    [mdlSFA] = sequential_lm(spikert,sideway,forward,angular);
    % reshape response and predictor data
    spikert = reshape(spikert,[],1);
    forward = reshape(forward,[],1);
    angular = reshape(angular,[],1);
    sideway = reshape(sideway,[],1);
    % store velocity for scatter plot
    allVelocity{nf,1} = forward;
    allVelocity{nf,2} = angular;
    allVelocity{nf,3} = sideway;
    allSpikert{nf,1} = spikert;

    % if sufficient data was present, continue model
    % else omit this fly
    predictor_t(1,nf) = t(2) * sum(~isnan(allVelocity{nf,1}));
    predictor_t(2,nf) = t(2) * sum(~isnan(allVelocity{nf,2}));
    predictor_t(3,nf) = t(2) * sum(~isnan(allVelocity{nf,3}));
    incFly(nf) = sum(predictor_t(:,nf)>cutoffD);

    if incFly(nf)==3
        % store R2 values
        fasR2(:,nf) = mdlFAS.R2;
        asfR2(:,nf) = mdlASF.R2;
        sfaR2(:,nf) = mdlSFA.R2;

        disp('Darkness behavior sufficient, continuing fit...')
        % second, fit visual input from quiescent data w/oscillating target
        thisFile = fullfile(pursuitFolder,join([thisFly '_int.mat'],""));
        % fetch pursuit dataset
        dfP = load(thisFile);

        % pull response data during QUIESCENT behavior only
        forward = dfP.int_forward(:,:,speedSelect);
        rest_idx = forward<restThreshold;
        spikert = dfP.int_spikert(:,:,speedSelect);
        rest_spikert = spikert;
        rest_spikert(~rest_idx) = nan;
        rest_spikert = reshape(rest_spikert,[],1);
        % pull predictor data and lag shift
        panelps = lagShift(dfP.int_panelps(:,:,speedSelect),t,obj_xcorr.Obj);
        % and determine panel direction and reshape
        paneldr = panelDirection(panelps);
        panelps = reshape(panelps,[],1);
        paneldr = reshape(paneldr,[],1);
        % normalize direction and pull IPSI only
        panelpsI = panelps;
        paneldrI = paneldr;
        paneldrI(panelpsI<0) = nan;
        panelpsI(panelpsI<0) = nan;

        % fit linear models for visual input
        predictorOrder2(1,:) = {'+Obj Position','+Obj Direction'};
        [mdlPD] = sequential_lm(rest_spikert,panelpsI,paneldrI,0);
        predictorOrder2(2,:) = {'+Obj Direction','+Obj Position'};
        [mdlDP] = sequential_lm(rest_spikert,paneldrI,panelpsI,0);
        % store R2 values
        pdR2(:,nf) = mdlPD.R2;
        dpR2(:,nf) = mdlDP.R2;
        % store obj position for scatter plot
        allPanelps{nf,1} = panelpsI;
        allPanelps{nf,2} = paneldrI;
        allSpikert{nf,2} = rest_spikert;
        % determine how much data was actually used for fit
        predictor_t(4,nf) = t(2) * sum(~isnan(allPanelps{nf,1}));

        % third, apply models to pursuit data
        % pull response data during PURSUIT behavior only
        [~,~,~,pur_spikert] = pursuitFinder(forward,0,0,spikert,dfP.int_time, pursThreshold);
        spikert = reshape(spikert,[],1);
        pur_spikert = reshape(pur_spikert,[],1);
        % pull predictor data and lag shift
        t = dfP.int_time;
        forward = lagShift(dfP.int_forward(:,:,1),t,settings.fwdLag);
        angular = lagShift(dfP.int_angular(:,:,1),t,settings.angLag);
        sideway = lagShift(dfP.int_sideway(:,:,1),t,settings.sidLag);
        % reshape
        forward = reshape(forward,[],1);
        angular = reshape(angular,[],1);
        sideway = reshape(sideway,[],1);
        % normalize direction and pull IPSI only
        angular(angular<0) = nan;
        sideway(sideway<0) = nan;
        % plot 2-3 trials depending on availability
        nTrial(nf,:) = size(dfP.int_spikert(:,:,speedSelect));
        eTime = dfP.int_time'; %time
        if nTrial(nf,2)<3
            pWin = 1:nTrial(1,1)*2;
            eTime = [eTime (eTime+eTime(end)+eTime(2))];
            eTime = reshape(eTime,[],1);
        else
            pWin = 1:nTrial(1,1)*3;
            eTime = [eTime (eTime+eTime(end)+eTime(2)) (eTime+eTime(end)*2+eTime(3))];
            eTime = reshape(eTime,[],1);
        end

        % compute predictions based on fits
        fwdFit = (mdlFAS.X1(1)*forward)+mdlFAS.Intercept(1);
        angFit = (mdlFAS.X1(2)*angular)+mdlFAS.Intercept(2);
        posFit = (mdlPD.X1(1)*panelpsI)+mdlPD.Intercept(1);
        dirFit = (mdlPD.X1(2)*paneldrI)+mdlPD.Intercept(2);
        % replace nans with 0s
        fwdFit(isnan(fwdFit)) = 0;
        angFit(isnan(angFit)) = 0;
        posFit(isnan(posFit)) = 0;
        dirFit(isnan(dirFit)) = 0;

        % assess fits
        predictorOrder(1,:) = {'+Forward','+Angular','+Obj Position','+Obj Direction'};
        fullR2_fapd(1,nf) = calculate_r_squared(pur_spikert,fwdFit);
        fullR2_fapd(2,nf) = calculate_r_squared(pur_spikert,fwdFit+angFit);
        fullR2_fapd(3,nf) = calculate_r_squared(pur_spikert,fwdFit+angFit+posFit);
        fullR2_fapd(4,nf) = calculate_r_squared(pur_spikert,fwdFit+angFit+posFit+dirFit);

        predictorOrder(2,:) = {'+Obj Position','+Obj Direction','+Forward','+Angular'};
        fullR2_pdfa(1,nf) = calculate_r_squared(pur_spikert,posFit);
        fullR2_pdfa(2,nf) = calculate_r_squared(pur_spikert,posFit+dirFit);
        fullR2_pdfa(3,nf) = calculate_r_squared(pur_spikert,posFit+dirFit+fwdFit);
        fullR2_pdfa(4,nf) = calculate_r_squared(pur_spikert,posFit+dirFit+fwdFit+angFit);

        % plot fits vs data
        figure; set(gcf,'Position',[50 50 1800 900]);
        tiledlayout(4,1,'TileSpacing','compact')
        sr_lim = [0 ceil(max(spikert))];

        nexttile; plot(eTime,panelps(pWin),'Color',p_colors{1}); axis tight
        yline(0); ylabel('Object Position (deg)'); legend('Object Pos')

        nexttile; hold on; plot(eTime,spikert(pWin),'k');
        plot(eTime,fwdFit(pWin),'Color',v_colors{1}); plot(eTime,angFit(pWin),'Color',v_colors{2}); axis tight
        ylim(sr_lim); legend('Firing Rate','Fwd Fit','Ang Fit')
        ylabel('Firing Rate')

        nexttile; hold on; plot(eTime,spikert(pWin),'k');
        plot(eTime,posFit(pWin),'Color',p_colors{1}); plot(eTime,dirFit(pWin),'Color',p_colors{2}); axis tight
        ylim(sr_lim); legend('Firing Rate','Target Pos Fit','Target Dir Fit')
        ylabel('Firing Rate')

        fullFit = posFit+dirFit+fwdFit+angFit;
        nexttile; hold on; plot(eTime,spikert(pWin),'k');
        plot(eTime,fullFit(pWin),'Color','#A2142F'); axis tight
        ylim(sr_lim); legend('Firing Rate','Full fit')
        ylabel('Firing Rate'); xlabel('Time (s)')

        % save image
        sgtitle([strrep(thisFly,'_',' ') ' model fit'])
        cd(plotFolder)
        plotname = join(['Fly' num2str(nf) '_modelfit'],"");
        saveas(gcf,[plotname '.png']);
    else
        disp('Darkness behavior insufficient, fly omitted.')
    end
end

%% plot data inclusion summary
figure; set(gcf,'Position',[50 50 500 900]);
tiledlayout(2, 1,'TileSpacing','compact')
nexttile
for nf = 1:nFliesD
    for v=1:3
        if incFly(nf)==3
            thisColor = v_colors{v};
        else
            thisColor = [1 1 1];
        end
        hold on; plot(nf,predictor_t(v,nf),'o','Color',v_colors{v},'MarkerFaceColor',thisColor)
    end
end
yline(25); xticks(1:nFliesD);
xlabel('Fly #'); ylabel('Time Available to Fit (s)')
legend(v_names)

nexttile
hold on; plot(1:nFliesD,predictor_t(4,:),'o','Color',p_colors{1},'MarkerFaceColor',p_colors{1})
yline(0); xticks(1:nFliesD);
xlabel('Fly #'); ylabel('Time Available to Fit (s)')
legend(p_names{1})
plotname = 'cutoff_predictors';
saveas(gcf,[plotname '.png']);


%% scatter velocity predictors and compare to model fit
% initialize
figure; set(gcf,'Position',[50 50 3500 900]);
tiledlayout(3, nFliesD,'TileSpacing','compact')
sr_lim = [0 160];

% for each velocity
for v = 1:3
    % for each fly
    for nf = 1:nFliesD
        % generate linear model
        dataTbl = table(allVelocity{nf,v},allSpikert{nf,1},'VariableNames',{v_names{v},'Spike Rate'});
        mdl = fitlm(dataTbl);

        nexttile
        f = plot(mdl);
        % get handles to plot components
        dataHandle = findobj(f,'DisplayName','Data');
        fitHandle = findobj(f,'DisplayName','Fit');
        cbHandles = findobj(f,'DisplayName','Confidence bounds');
        cbHandles = findobj(f,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
        % adjust colors
        if incFly(nf)==3
            this_color = v_colors{v};
            this_trial = trial_grey;
        else
            this_color = trial_grey - 0.1;
            this_trial = trial_grey + 0.2;
        end
        fitHandle.Color = this_color; fitHandle.LineWidth = 3; set(cbHandles,'Color','k')
        dataHandle.Color = this_trial; dataHandle.Marker = '.'; dataHandle.MarkerSize = 1;
        axis tight; ylim(sr_lim); title(['Rsqd = ' num2str(mdl.Rsquared.Adjusted)]);
        if nf>1
            ylabel('')
        end
    end
end

% save
cd(plotFolder)
plotname = 'scatter_predictors_vel';
saveas(gcf,[plotname '.png']);

%% scatter visual target predictors and compare to model fit
% initialize
figure; set(gcf,'Position',[50 50 3500 600]);
tiledlayout(2, nFliesD,'TileSpacing','compact')
p_range = [0 40; -1.5 1.5];
sr_lim = [0 150];

% for each target feature
for p = 1:2
    % for each fly
    for nf = 1:nFliesD
        nexttile
        if incFly(nf)==3
            % generate linear model
            dataTbl = table(allPanelps{nf,p},allSpikert{nf,2},'VariableNames',{p_names{p},'SpikeRate'});
            mdl = fitlm(dataTbl);
            this_color = p_colors{p};

            f = plot(mdl);
            % get handles to plot components
            dataHandle = findobj(f,'DisplayName','Data');
            fitHandle = findobj(f,'DisplayName','Fit');
            cbHandles = findobj(f,'DisplayName','Confidence bounds');
            cbHandles = findobj(f,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
            % adjust colors
            fitHandle.Color = this_color; fitHandle.LineWidth = 3; set(cbHandles,'Color','k')
            dataHandle.Color = trial_grey; dataHandle.Marker = '.'; dataHandle.MarkerSize = 1;
            axis tight; ylim(sr_lim); xlim(p_range(p,:))
            title(['Rsqd = ' num2str(mdl.Rsquared.Adjusted)]);
        end
    end
end
% save
cd(plotFolder)
plotname = 'scatter_predictors_obj';
saveas(gcf,[plotname '.png']);

%% summarize R2 for velocity predictors
figure; set(gcf,'Position',[100 100 1200 400]);
tiledlayout(1,4,'TileSpacing','compact')
r2_lim = [0 0.8];
% for each run, plot individual R2
nRun = 3;
for m = 1:nRun
    % select stored R2 values
    switch m
        case 1
            r2 = fasR2;
        case 2
            r2 = asfR2;
        case 3
            r2 = sfaR2;
    end
    % calculate mean, SEM
    r2_mean = median(r2,2,'omitnan');
    r2_sem = std(r2,[],2,'omitnan')./sqrt(sum(incFly==3));

    nexttile
    hold on; plot(1:nRun,r2,':o','Color',trial_grey)
    errorbar(1:nRun,r2_mean,r2_sem,'-o','Color','k')
    axis padded; ylim(r2_lim)
    ylabel('R2'); xlabel('Predictor')
    xticks(1:nRun); xticklabels(predictorOrder1(m,:)); xtickangle(45)
end

% generate interaction comparison
interactR2 = nan(2,nFliesD);
for nf = 1:nFliesD
    if incFly(nf)==3
        % fit forward+angular
        dataTbl = table(allVelocity{nf,1:2},allSpikert{nf,1},'VariableNames',{v_names{1:2},'Spike Rate'});
        mdl1 = fitlm(dataTbl,'interactions');
        % fit forward+sideways
        dataTbl = table(allVelocity{nf,[1 3]},allSpikert{nf,1},'VariableNames',{v_names{[1 3]},'Spike Rate'});
        mdl2 = fitlm(dataTbl,'interactions');
        % store
        interactR2(1,nf) = mdl1.Rsquared.Adjusted;
        interactR2(2,nf) = mdl2.Rsquared.Adjusted;
    end
end
% plot
interactMean = median(interactR2,2,'omitnan');
interactSEM = std(interactR2,[],2,'omitnan')./sqrt(sum(incFly==3));
nexttile
hold on; plot(1:2,interactR2,'-o','Color',trial_grey)
errorbar(1:2,interactMean,interactSEM,'-o','Color','k')
axis padded; ylim(r2_lim)
ylabel('R2'); xlabel('Predictor Interaction')
xticks(1:2); xticklabels({'Fwd+Ang','Fwd+Sid'}); xtickangle(45)

% save image
sgtitle('fit to velocity when fly walking in darkness w P1')
cd(plotFolder)
plotname = 'R2_predictors_vel';
saveas(gcf,[plotname '.png']);

%% summarize R2 for visual target predictors
figure; set(gcf,'Position',[100 100 1000 400]);
tiledlayout(1,3,'TileSpacing','compact')
r2_lim = [0 1];
% for each run
nRun = 2;
for m = 1:nRun
    % select stored R2 values
    switch m
        case 1
            r2 = pdR2;
        case 2
            r2 = dpR2;
    end
    % calculate mean, SEM
    r2_mean = median(r2,2,'omitnan');
    r2_sem = std(r2,[],2,'omitnan')./sqrt(sum(incFly==3));

    nexttile
    hold on; plot(1:nRun,r2,':o','Color',trial_grey)
    errorbar(1:nRun,r2_mean,r2_sem,'-o','Color','k')
    axis padded; ylim(r2_lim)
    ylabel('R2'); xlabel('Predictor')
    xticks(1:nRun); xticklabels(predictorOrder2(m,:)); xtickangle(45)
end

% generate interaction comparison
interactR2 = nan(1,nFliesD);
for nf = 1:nFliesD
    if incFly(nf)==3
        % fit forward+angular
        dataTbl = table(allPanelps{nf,1:2},allSpikert{nf,2},'VariableNames',{p_names{1:2},'Spike Rate'});
        mdl1 = fitlm(dataTbl,'interactions');
        % store
        interactR2(1,nf) = mdl1.Rsquared.Adjusted;
    end
end
% plot
interactMean = median(interactR2,2,'omitnan');
interactSEM = std(interactR2,[],2,'omitnan')./sqrt(sum(incFly==3));
nexttile
hold on; plot(1,interactR2,'-o','Color',trial_grey)
errorbar(1,interactMean,interactSEM,'-o','Color','k')
axis padded; ylim(r2_lim)
ylabel('R2'); xlabel('Predictor Interaction')
xticks(1); xticklabels('Pos+Dir'); xtickangle(45)

sgtitle('fit to visual input when fly at rest during pursuit trials')
% save image
cd(plotFolder)
plotname = 'R2_predictors_obj';
saveas(gcf,[plotname '.png']);

%% plot model performance
figure; set(gcf,'Position',[100 100 800 400]);
% for each run
nRun = 2;
r2_lim = [-1 1];
for m = 1:nRun
    % select stored R2 values
    switch m
        case 1
            r2 = fullR2_fapd;
        case 2
            r2 = fullR2_pdfa;
    end
    % calculate mean, SEM
    r2_mean = mean(r2,2,'omitnan');
    r2_sem = std(r2,[],2,'omitnan')./sqrt(nFliesD);

    subplot(1,nRun,m)
    hold on; plot(1:4,r2,':o','Color',trial_grey)
    errorbar(1:4,r2_mean,r2_sem,'-o','Color','k')
    axis padded; ylim(r2_lim); yline(0)
    ylabel('R2'); xlabel('Predictor')
    xticks(1:4); xticklabels(predictorOrder(m,:)); xtickangle(45)
end
sgtitle('model performance on pursuit data')
% save image
cd(plotFolder)
plotname = 'R2_FULL';
saveas(gcf,[plotname '.png']);
