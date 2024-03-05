%controlsimulation
% model behavior of a PD controller when presented with an oscillating
% target
%
% INPUT
% sweepSpeed        - speed of target (e.g., 50 70 100deg/s)
% controllerDelay   - controller delay (e.g. 300msec.)
% stimType          - 1 for sawtooth, 2 for pseduorandom
%                       note: can be removed, but some problems with lag
% plotSim           - 1 to plot simulation, 0 to not
%
% OUTPUT
% cc_pd             - correlation coefficient for model performance
% amp               - amplitude of each output at peak
% phase             - delay of each output relative to target motion
%
% CREATED - 09/12/2022 MC and RW
% UPDATED - 09/19/2022 MC removed weight variable, adjusted derivative amp
%           09/20/2022 MC added delay calculation
%           02/07/2023 MC added delay output
%

function [cc_pd,amp,phase] = controlsimulation(stimwave,funcFreq,controllerDelay,stimType,plotSim)
%% initialize
disp('Starting control simulation...')
% set time array
duration = length(stimwave)/funcFreq;
t = 1/funcFreq:1/funcFreq:duration;


%% calculate error and output for each controller, then combine

% smooth stimwave input
sm_stimwave = smoothdata(stimwave,'gaussian',funcFreq);
sm_stimwave = (sm_stimwave./range(sm_stimwave./2)).*range(stimwave./2); %normalize
%sm_stimwave=stimwave;
% generate delay array
delay = zeros(1,round((controllerDelay/1000)*funcFreq));

% generate position unit
pEC = sm_stimwave;                          %position error = obj position
% add delay
pEC = [delay pEC];
pEC = pEC(1:length(t));
% normalize position error unit
pEC = range(stimwave)*pEC/range(pEC);

% generate derivative error unit
dEC = diff([sm_stimwave nan])/(1/funcFreq);  %derivative error = derivative/change t
% add delay
dEC = [delay dEC];
dEC = dEC(1:length(t));
% normalize
dEC = dEC .* 0.35;

% generate combined PD error unit
cEC = pEC + dEC;
% calculate correlation coefficient between waveform input and PD output
cc_pd = corrcoef(stimwave,cEC);


%% find peaks to determine output phase relative to stimulus

% initialize
phase={};

if stimType==1
    pdis = 300; %min peak distance
    phght = 1; %min peak height
    [~,idx_s] = findpeaks(stimwave,'MinPeakHeight',phght,'MinPeakDistance',pdis); %stimulus
    [pk_p,idx_p] = findpeaks(pEC,'MinPeakHeight',phght,'MinPeakDistance',pdis); %pos
    [pk_d,idx_d] = findpeaks(dEC,'MinPeakHeight',phght,'MinPeakDistance',pdis); %der
    [pk_c,idx_c] = findpeaks(cEC,'MinPeakHeight',phght,'MinPeakDistance',pdis);

    if length(idx_p)>1
        % remove first/last peaks
        idx_s = idx_s(2:end);
        idx_p = idx_p(2:end);
        idx_d = idx_d(2:end-1);
        idx_c = idx_c(2:end-1);
        % calculate mean phase lag between each output and the stimulus
        phase.pos = mean((t(idx_p) - t(idx_s)))*1000;
        phase.der = mean((t(idx_d) - t(idx_s)))*1000;
        phase.combo = mean((t(idx_c) - t(idx_s)))*1000;
        % calculate mean output magnitude for each output
        amp.pos = mean(pk_p(2:end));
        amp.der = mean(pk_d(2:end-1));
        amp.combo = mean(pk_c(2:end-1));
    end
else
    idx_p = [];
end


%% find offset

cEC_n = (cEC./range(cEC./2)).*range(stimwave./2); %normalize
offset = sm_stimwave-cEC_n;

%% plot stimulus and controller

if plotSim
    % generate figure
    figure(1); clf;
    if duration>15
        fwidth = 1500; %max figure size
    else
        fwidth = 100*duration;
    end
    set(gcf,'Position',[100 100 fwidth 800])
    lw = 1;

    
    s1 = subplot(2,1,1);
    % plot target waveform
    plot(t,stimwave,'Color','k','LineWidth',lw)
    yline(0,'Color','#7E2F8E')
    % plot output peaks
    if ~isempty(idx_p)
        xline(t(idx_s),'Color','k')
        xline(t(idx_p),'Color','#D95319')
        xline(t(idx_d),'Color','#77AC30')
        xline(t(idx_c),'Color','#0072BD')

        dim = [0.8 0.6 0.3 0.3];
        str = {['P(delay) = ' num2str(phase.pos)],['D(delay) = ' num2str(phase.der)],['PD(delay) = ' num2str(phase.combo)]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
    ylabel('Target Position (deg)');
    axis tight


    s2 = subplot(2,1,2);
    hold on
    % plot position and derivative controller outputs
    plot(t,pEC,'Color','#D95319','LineWidth',lw)
    plot(t,dEC,'Color','#77AC30','LineWidth',lw)
    % plot combined PD controller output
    plot(t,cEC,'Color','#0072BD','LineWidth',lw*2)
    yline(0,'Color','k')
    hold off

    axis tight
    ylabel('Controller Output (deg/s)');
    
    
%     s3 = subplot(3,1,3);
%     % plot offset
%     area(t,offset,'FaceColor','#0072BD','FaceAlpha',.3,'EdgeAlpha',0)
%     axis tight
%     ylabel('Offset (deg))');

    linkaxes([s1,s2],'x');
    xlabel('Time (s)')

end

disp('control simulation complete.')
