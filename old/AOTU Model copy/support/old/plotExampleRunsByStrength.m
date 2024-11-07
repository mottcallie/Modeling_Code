function plotExampleRunsByStrength(kValues, runSettings, predicted_RF, noiseLevel, startPos, simDuration, thisSynapse, alpha, shift)
    % plotExampleRunsByStrength
    % This function plots one example run from the strength 0 and strength 1 condition for each k.
    % kValues: Array of steering gain (k) values
    % runSettings: Settings structure for the model run
    % predicted_RF: Tuning structure for the model
    % noiseLevel: Noise level for the simulation
    % startPos: Starting position for the simulation
    % simDuration: Duration of the simulation
    % thisSynapse: Synapse type (e.g., inhibitory)
    % alpha: Alpha value for the current simulation
    % shift: Shift value for the current simulation

    % Define strength values for comparison (Strength 1 and Strength 0)
    strengthValues = [1, 0];
    nK = length(kValues);  % Number of k values

    % Loop over k values (steering gain)
    figure;
    tiledlayout(nK, 2, 'TileSpacing', 'compact');  % Plot one row for each k, two columns (strength 1 and strength 0)

    for kIdx = 1:nK
        thisK = kValues(kIdx);

        for strengthIdx = 1:2
            thisStrength = strengthValues(strengthIdx);

            % Adjust tuning for the current strength
            thisTuning = predicted_RF;
            thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;  % Adjust tuning by strength

            % Run the AOTU steering model with adjusted tuning
            [timebase, visobj_history, input_history, rotvel_history] = ...
                aotu_steering_model(thisTuning, noiseLevel, startPos, thisSynapse, simDuration, thisK, alpha, runSettings);

            % Plot results for the current strength and k
            nexttile;
            plot(timebase, visobj_history);
            xlabel('Time (s)');
            ylabel('Object Position (deg)');
            if strengthIdx == 1
                title(['Strength 1, k = ' num2str(thisK)]);
            else
                title(['Strength 0, k = ' num2str(thisK)]);
            end
        end
    end

    % Save the figure
    sgtitle(['Example Runs for Strength 1 vs Strength 0 Across k']);
    saveas(gcf, fullfile(runSettings.folder.summary, 'Example_Run_Strength1_vs_0.png'));
end
