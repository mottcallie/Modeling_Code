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
%
% CREATED - 09/12/2022 MC and RW
% UPDATED - 09/19/2022 MC removed weight variable, adjusted derivative amp
%           09/20/2022 MC added delay calculation
%

function [cc_pd] = controlsimulation(stimwave,funcFreq,controllerDelay,stimType,plotSim)
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

if stimType==1
    pdis = 300; %min peak distance
    phght = 1; %min peak height
    [~,is] = findpeaks(stimwave,'MinPeakHeight',phght,'MinPeakDistance',pdis); %stimulus
    [~,ip] = findpeaks(pEC,'MinPeakHeight',phght,'MinPeakDistance',pdis); %pos
    [~,id] = findpeaks(dEC,'MinPeakHeight',phght,'MinPeakDistance',pdis); %der
    [~,ic] = findpeaks(cEC,'MinPeakHeight',phght,'MinPeakDistance',pdis);

    if length(ip)>1
        % remove first/last, as these outputs are often distorted
        is(1) = [];
        ip(1) = [];
        id([1 end])=[];
        ic([1 end])=[];

        pDelay = mean((t(ip) - t(is)))*1000;
        dDelay = mean((t(id) - t(is)))*1000;
        cDelay = mean((t(ic) - t(is)))*1000;
    end
else
    ip = [];
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
    if ~isempty(ip)
        xline(t(is),'Color','k')
        xline(t(ip),'Color','#D95319')
        xline(t(id),'Color','#77AC30')
        xline(t(ic),'Color','#0072BD')

        dim = [0.8 0.6 0.3 0.3];
        str = {['P(delay) = ' num2str(pDelay)],['D(delay) = ' num2str(dDelay)],['PD(delay) = ' num2str(cDelay)]};
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
