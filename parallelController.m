% parallelController
% modeling code for building and testing a simple input-output controller
% composed of 2 systems in parallel. The first is a "generalized" system
% with activity that is proportional to current error with low gain output
% drive. The second is a "dedicated" system with activity that is driven by
% error near the midline only with high gain output drive.
%
%% initialize

% generate system receptive fields
% region of visual space covered by generalized system
rf_general = 90; %deg
% region of visual space covered by specialized system
rf_dedicated  = 45; %deg

% region of binocular overlap
binocOverlapSize = 15; %deg
% direction selectivity (ipsiversive>contraversive)
dsIndex = 1.2;

%% plot input-output curve for both systems

% initialize
figure;
sampleVisualInput = -90:90;

[out_gen] = generalSystem(sampleVisualInput,rf_general,1);
subplot(1,3,1)
plot(sampleVisualInput,out_gen)
xline(0);yline(0);

[out_ded] = dedicatedSystem(sampleVisualInput,rf_dedicated,0.8);
subplot(1,3,2)
plot(sampleVisualInput,out_ded)
xline(0);yline(0);

subplot(1,3,3)
out_combo = smooth(out_gen+out_ded,100);
plot(sampleVisualInput,out_combo)
xline(0);yline(0);
