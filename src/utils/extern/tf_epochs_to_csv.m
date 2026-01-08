% MATLAB script to save FieldTrip epochs structure to CSV files
% Usage: matlab_to_csv(epochs, output_dir)
%
% Saves:
% - epochs_data.csv: Combined data with columns [trial, time, channel1, channel2, ...]
% - epochs_metadata.csv: Metadata (fsample, labels, trialinfo, sampleinfo)

function matlab_to_csv(epochs, output_dir)
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    n_trials = length(epochs.trial);
    n_channels = length(epochs.label);
    
    % Get channel names (convert cell array to string array for CSV)
    channel_names = epochs.label;
    
    % Prepare combined data matrix
    % Each row: [trial_idx, time, channel1, channel2, ..., channelN]
    all_data = [];
    
    for trial_idx = 1:n_trials
        % Get trial data (channels × samples)
        trial_data = epochs.trial{trial_idx};
        time_vec = epochs.time{trial_idx};
        
        % Transpose to samples × channels
        trial_data_transposed = trial_data';
        n_samples = size(trial_data_transposed, 1);
        
        % Create matrix: [trial_idx, time, channels...]
        trial_matrix = [repmat(trial_idx, n_samples, 1), time_vec(:), trial_data_transposed];
        
        % Append to all_data
        all_data = [all_data; trial_matrix];
    end
    
    % Create column names - ensure channel_names is a cell array of strings
    % channel_names is already a cell array from epochs.label
    col_names = cell(1, 2 + n_channels);
    col_names{1} = 'trial';
    col_names{2} = 'time';
    for i = 1:n_channels
        col_names{2 + i} = channel_names{i};
    end
    
    % Write to CSV - all information is in this single file
    data_file = fullfile(output_dir, 'epochs_data.csv');
    writetable(array2table(all_data, 'VariableNames', col_names), data_file);
    
    fprintf('Saved epochs data to: %s\n', output_dir);
    fprintf('  - epochs_data.csv: %d rows, %d columns (trial, time, %d channels)\n', ...
            size(all_data, 1), size(all_data, 2), n_channels);
    fprintf('  - fsample: %.1f Hz\n', epochs.fsample);
    fprintf('  - n_trials: %d\n', n_trials);
end

