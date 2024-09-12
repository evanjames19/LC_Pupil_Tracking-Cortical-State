% Define main folder and file names
root_folder = 'C:\Users\Neuron\Documents\MATLAB\EEG_Analysis\LF_LC_stim_analysis\LC_PAPER\ANIMALS_STIM';
file_name = 'stim_onsets.mat';
pupil_name = 'pupilTrace.mat';
trial_t = 'trial_type.mat';
v19_name = 'low-A-019.dat';
warning('off');

% List all animal ID folders
animal_folders = dir(root_folder);
animal_folders = animal_folders([animal_folders.isdir]);
animal_folders = animal_folders(~ismember({animal_folders.name}, {'.', '..'}));

% Sampling frequency and frequency bands
Fs = 1000;
freq_bands = [1 5; 5 10; 8 12; 13 30; 30 100];

% Pupil dilation detection parameters
windowSize = 1000;
zscoreThreshold = 1;
maxDilations = 1000;

% Main structure to store all session data
main_struct = [];
%%
for i = 1:length(animal_folders)
    animal_id = animal_folders(i).name;
    animal_path = fullfile(root_folder, animal_id);
    fprintf('Processing animal folder: %s\n', animal_path); % Debugging: Confirm animal folder path
    
    % List all session date folders for the current animal
    session_folders = dir(animal_path);
    session_folders = session_folders([session_folders.isdir]);
    session_folders = session_folders(~ismember({session_folders.name}, {'.', '..'}));
    
    %%
    % Loop through each session folder
    for q = 1:length(session_folders)
        current_folder = fullfile(animal_path, session_folders(q).name);
        fprintf('Processing session folder: %s\n', current_folder); % Debugging: Confirm session folder path
        
        % Define file paths
        pupil_file = fullfile(current_folder, 'pupilTrace.mat');
        current_file = fullfile(current_folder, 'stim_onsets.mat');
        v19_file = fullfile(current_folder, 'low-A-019.dat');
        trial_file = fullfile(current_folder, 'trial_type.mat');

        animal_double = str2double(animal_id);
        animal = animal_double;
        if isfile(pupil_file) && isfile(trial_file) && isfile(v19_file) && isfile(current_file)
            fprintf('Processing pupil file: %s\n', pupil_file);
            fprintf('Processing EEG file: %s\n', v19_file);
            fprintf('Processing trial file: %s\n', trial_file);
            fprintf('Processing stim onset file: %s\n', current_file);
            stimulation_times = load(current_file);
            stim_times = stimulation_times.stimOnset;
            trialstypes = load(trial_file);
            trial_types = trialstypes.trial_type;
            
            % Load the data
            pupilData = load(pupil_file);
            trialData = load(trial_file);

            fileinfo = dir(v19_file);
            num_samples = fileinfo.bytes / 2; % int16 = 2 bytes
            fid = fopen(v19_file, 'r');
            v19 = fread(fid, num_samples, 'int16');
            fclose(fid);
            v19 = v19 * 0.195; % convert to microvolts
            eegData = v19(1001:end); % Exclude the first 1000 samples if needed
            
            eegData = zscore(eegData);
            eegData = lowpass(eegData, 100, 1000);
            
            % Extract pupil trace and calculate z-scores
            pupilTrace = pupilData.pupilTrace;
            pupilTraceZScores = zscore(pupilTrace);
            pupil_diff = [diff(pupilTraceZScores); 0];

            % Time vector for original data (10 Hz)
            originalTime = linspace(0, length(pupilTraceZScores) / 10, length(pupilTraceZScores));

            % New time vector for 100 Hz
            newTime = linspace(0, length(pupilTraceZScores) / 10, length(pupilTraceZScores) * 10);

            % Spline Interpolation
            interpolatedData = interp1(originalTime, pupilTraceZScores, newTime, 'spline');
            interplotated_diff = interp1(originalTime, pupil_diff, newTime, 'spline');

            pupilData = interpolatedData';
            differential_pupil_data = interplotated_diff';

            Fs = 10; % Sampling frequency
            cutoff_frequency = 0.1; % Cutoff frequency for the high-pass filter
            
            % Design a 2rd-order high-pass Butterworth filter
            [B, A] = butter(2, cutoff_frequency / (Fs / 2), 'high');
            
            % Apply the filter to your data using filtfilt for zero-phase filtering
            pupilTraceZScores1 = filtfilt(B, A, pupilData);

            % Initialize the pupil structure for this session
            pupil_struct = repmat(struct('Date', [], 'Animal', [], 'RawPupil', [], 'PupilIndexes', [], ...
                'StimulationFrequency', [], 'DeltaPwr', [], 'ThetaPwr', [], 'AlphaPwr', [], 'BetaPwr', [], ...
                'LowGammaPwr', [], 'HighGammaPwr', [], 'Dilation', [], 'PreDilationSpectrum', [], ...
                'DuringDilationSpectrum', [], 'PostDilationSpectrum', []), maxDilations, 1);

            date = session_folders(q).name;
            date = str2double(date);

            structIndex = 1;
            i = 1;

            baselineStabilityThreshold = 200; % Number of points for baseline stability
            peakReturnThreshold = 600; % Maximum number of points to return to baseline after peak
            zscoreThreshold = 1; % Define the threshold for z-score 
            % (e.g., 2 standard deviations)
            maxDilations = 1000; % Maximum number of dilations to detect
            detectedDilations = 0; % Initialize counter for detected dilations
            structIndex = 1; % Index for storing dilation events
            windowSize = 1000; % Size of the window to check for dilations

            %%
            while i <= length(pupilData) - windowSize + 1
                window = pupilData(i:i + windowSize - 1);

                % Check for baseline stability
                if all(abs(window(1:baselineStabilityThreshold)) < zscoreThreshold)
                    % Check for a significant jump in z-score
                    peakIndices = find(window > zscoreThreshold);
                    if ~isempty(peakIndices)
                        firstPeakIndex = peakIndices(1);

                        % Check if the peak is within the first part of the window to ensure 100 points
                        if firstPeakIndex <= windowSize - peakReturnThreshold
                            % Check for the return to baseline
                            postPeakIndices = find(window(firstPeakIndex:firstPeakIndex + peakReturnThreshold - 1) < zscoreThreshold);
                            if ~isempty(postPeakIndices)
                                % Adjust dilationEnd to ensure it's 100 points from the start
                                dilationEnd = i + windowSize - 1;

                                % Calculate the EEG data indices for the current dilation window
                                stim_start = i - 1; % Ensure start is within bounds
                                stim_end = dilationEnd; % Ensure end is within bounds

                                if (stim_start > 0)
                                    fprintf('Extracting EEG data from index %d to %d\n', stim_start * 10, stim_end * 10);
                                    stim_window = eegData(stim_start * 10:stim_end * 10 - 1);
                                    % Assign other dilation-related data to the structure
                                    pupil_struct(structIndex).Date = date;
                                    pupil_struct(structIndex).Animal = animal;
                                    pupil_struct(structIndex).StimulationFrequency = 0;
                                    fprintf('Extracting Pupil data from index %d to %d\n', i, dilationEnd);

                                    pupil_struct(structIndex).PupilIndexes = [i - 1, dilationEnd];
                                    pupil_struct(structIndex).RawPupil = pupilData(i:dilationEnd);
                                    
                                    % Define window size in samples
                                    windowSize = 1000; % 1 second window at 1000 Hz

                                    % Total number of windows
                                    totalWindows = 1000; % Set to 1000 windows

                                    % Initialize the output matrix
                                    windowMatrix = zeros(windowSize, totalWindows);

                                    % Calculate step size based on desired number of windows
                                    stepSize = (length(stim_window) - windowSize) / (totalWindows - 1);

                                    % Extract windows
                                    for w = 1:totalWindows
                                        startIndex = round((w - 1) * stepSize) + 1;
                                        endIndex = startIndex + windowSize - 1;
                                        windowMatrix(:, w) = stim_window(startIndex:endIndex);
                                    end

                                    % Frequency bands
                                    freq_bands = [1 3; 4 7; 8 12; 13 30; 30 55; 65 99];
                                    num_bands = size(freq_bands, 1);

                                    % Initialize matrix to store average band power (5 frequency bands x 1000 windows)
                                    avg_band_power = zeros(num_bands, size(windowMatrix, 2));

                                    % Sampling frequency
                                    fs = 1000; % 1000 Hz

                                    % Frequency vector
                                    N = size(windowMatrix, 1);
                                    f = (0:(N/2)-1) * (fs/N);

                                    % Initialize PSD matrix
                                    psd_matrix = zeros(size(windowMatrix, 2), length(f));

                                    % Loop through each window
                                    for k = 1:size(windowMatrix, 2)
                                        % Compute the PSD for the current window
                                        current_window = windowMatrix(:, k);
                                        window_psd = abs(fft(current_window)).^2 / N;
                                        psd_vals = window_psd(1:N/2); % Only consider positive frequencies
                                        psd_matrix(k, :) = psd_vals; % Store the PSD

                                        % Loop through each frequency band
                                        for p = 1:num_bands
                                            % Extract the indices corresponding to the current frequency band
                                            band_indices = f >= freq_bands(p, 1) & f <= freq_bands(p, 2);

                                            % Compute the average power in the current frequency band
                                            avg_band_power(p, k) = trapz(f(band_indices), psd_vals(band_indices));
                                        end
                                    end

                                    % Assign EEG Pwr to struct;
                                    pupil_struct(structIndex).DeltaPwr = avg_band_power(1, :)';
                                    pupil_struct(structIndex).ThetaPwr = avg_band_power(2, :)';
                                    pupil_struct(structIndex).AlphaPwr = avg_band_power(3, :)';
                                    pupil_struct(structIndex).BetaPwr = avg_band_power(4, :)';
                                    pupil_struct(structIndex).LowGammaPwr = avg_band_power(5, :)';
                                    pupil_struct(structIndex).HighGammaPwr = avg_band_power(6, :)';
                                    pupil_struct(structIndex).Dilation = 1;
                                    pupil_struct(structIndex).PSD = psd_matrix; % Save the PSD

                                    % Increment the structIndex to move to the next dilation event
                                    fprintf('Processed dilation #%d, structIndex: %d\n', structIndex, structIndex); % Debug: Print after processing each dilation

                                    structIndex = structIndex + 1;
                                    i = dilationEnd; % Move to the end of the current dilation

                                    continue;
                                end
                            end
                        end
                    end
                end
                i = i + 1; % Move to the next point if no dilation is found
            end

            temp_new_date_struct = [];

            stimulation_onset = stim_times; % Your stimulation onset vector
            num_stimulations = length(stimulation_onset);
            extracted_windows = cell(1, num_stimulations);
            sampling_rate = 100; % Pupil trace sampling rate (in fps)
            window_duration = 12; % Window duration in seconds

            % Get unique dates as integers
            spontaneousDates = [pupil_struct.Date];
            uniqueDates = unique([spontaneousDates]);

            current_file_dates = pupil_struct([pupil_struct.Date] == date);
            current_date_trials_size = length(current_file_dates);

            time_v = linspace(0, length(pupilData) / 100, length(pupilData));
            time_vector = time_v; % Your time vector data

            all_pupiltrace_indexes = pupil_struct.PupilIndexes;

            for k = 1:current_date_trials_size

                pupil_index = current_file_dates(k).PupilIndexes;
                pupil_vector = current_file_dates(k).RawPupil;
                proper_index_vector = pupil_index(1):pupil_index(2) - 1;

                adjust_onsets = stimulation_onset * 100;
                diff1 = abs(adjust_onsets - pupil_index(1));

                % Find the index of the minimum difference
                [~, minIndex] = min(diff1);
                closestValue = stimulation_onset(minIndex) * 100;

                if(abs(closestValue - pupil_index(1)) <= 1000)
                    [getMaxPupilTrace, idxMaxPupilTrace] = max(current_file_dates(k).RawPupil); % max value
                    extract_proper_val = proper_index_vector(idxMaxPupilTrace);

                    if((extract_proper_val - closestValue) >= 0)
                        difference = (extract_proper_val - closestValue);
                        fprintf('Identified dilation #%d overlap\n', k); % Debug: Print after processing each dilation
                        fprintf('Difference between proper and closest: #%d \n ', difference);
                        new_pupil_vector = pupilData(closestValue - 1000:(closestValue + 999));
                        current_file_dates(k).PupilIndexes = [round(closestValue) - 1000, round(closestValue) + 999];
                        current_file_dates(k).RawPupil = new_pupil_vector;
                        current_freq = trial_types(minIndex);
                        current_file_dates(k).StimulationFrequency = current_freq;

                        stim_start = closestValue * 10 - 10000;
                        stim_end = closestValue * 10 + 9999;

                        if(stim_end < length(eegData))
                            stim_window = eegData(stim_start:stim_end);

                            % Define window size in samples
                            windowSize = 1000; % 1 second window at 1000 Hz

                            % Total number of windows
                            totalWindows = 2000; % Set to 1000 windows

                            % Initialize the output matrix
                            windowMatrix = zeros(windowSize, totalWindows);

                            % Calculate step size based on desired number of windows
                            stepSize = (length(stim_window) - windowSize) / (totalWindows - 1);

                            % Extract windows
                            for w = 1:totalWindows
                                startIndex = round((w - 1) * stepSize) + 1;
                                endIndex = startIndex + windowSize - 1;
                                windowMatrix(:, w) = stim_window(startIndex:endIndex);
                            end

                            % Frequency bands
                            freq_bands = [1 3; 4 7; 8 12; 13 30; 30 55; 65 99];
                            num_bands = size(freq_bands, 1);

                            % Initialize matrix to store average band power (5 frequency bands x 1000 windows)
                            avg_band_power = zeros(num_bands, size(windowMatrix, 2));

                            % Sampling frequency
                            fs = 1000; % 1000 Hz

                            % Frequency vector
                            N = size(windowMatrix, 1);
                            f = (0:(N/2)-1) * (fs/N);

                            % Initialize PSD matrix
                            psd_matrix = zeros(size(windowMatrix, 2), length(f));

                            % Loop through each window
                            for b = 1:size(windowMatrix, 2)
                                % Compute the PSD for the current window
                                current_window = windowMatrix(:, b);
                                window_psd = abs(fft(current_window)).^2 / N;
                                psd_vals = window_psd(1:N/2); % Only consider positive frequencies
                                psd_matrix(b, :) = psd_vals; % Store the PSD

                                % Loop through each frequency band
                                for p = 1:num_bands
                                    % Extract the indices corresponding to the current frequency band
                                    band_indices = f >= freq_bands(p, 1) & f <= freq_bands(p, 2);

                                    % Compute the average power in the current frequency band
                                    avg_band_power(p, b) = trapz(f(band_indices), psd_vals(band_indices));
                                end
                            end

                            % Assign EEG Pwr to struct;
                            current_file_dates(k).DeltaPwr = avg_band_power(1, :)';
                            current_file_dates(k).ThetaPwr = avg_band_power(2, :)';
                            current_file_dates(k).AlphaPwr = avg_band_power(3, :)';
                            current_file_dates(k).BetaPwr = avg_band_power(4, :)';
                            current_file_dates(k).LowGammaPwr = avg_band_power(5, :)';
                            current_file_dates(k).HighGammaPwr = avg_band_power(6, :)';
                            current_file_dates(k).Dilation = 1;
                            current_file_dates(k).PSD = psd_matrix; % Save the PSD
                        end
                    end
                end
            end

            temp_new_date_struct = [temp_new_date_struct; current_file_dates];
            main_struct = [main_struct; temp_new_date_struct(1:structIndex-1)]; % Concatenate the structures
        else
            fprintf('One or more files not found in: %s\n', current_folder);
        end
    end
    spontaneous_struct = main_struct;
end
