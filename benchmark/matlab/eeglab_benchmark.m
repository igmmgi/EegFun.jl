function eeglab_benchmark(data_dir, n_files_to_process, run_ica_flag)
    % EEGLAB benchmark script for ERP pipeline
    % Ensure EEGLAB is in your MATLAB path before running.
    
    if nargin < 1
        error('Please provide the data directory as an argument (e.g., eeglab_benchmark(''/path/to/data''))');
    end
    
    if nargin < 2
        n_files_to_process = 0;
    end
    
    if nargin < 3
        run_ica_flag = true;
    else
        if ischar(run_ica_flag) || isstring(run_ica_flag)
            run_ica_flag = strcmpi(run_ica_flag, 'true') || strcmpi(run_ica_flag, '1');
        end
    end

    fprintf('Benchmarking EEGLAB pipeline...\n');
    t_start = tic;

    % Initialize EEGLAB paths and headless mode
    eeglab nogui;
    
    % Get the absolute path of the directory containing this script
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    
    if ~exist(data_dir, 'dir')
        error('Data directory not found: %s', data_dir);
    end
    file_list = dir(fullfile(data_dir, '*.bdf'));
    test_files = fullfile({file_list.folder}, {file_list.name});
    test_files = sort(test_files);
    
    if n_files_to_process > 0
        test_files = test_files(1:min(length(test_files), n_files_to_process));
    end

    for i = 1:length(test_files)
        % 1. Read
        EEG = pop_biosig(test_files{i});

        % Load the pre-compiled native EEGLAB 3D coordinate file
        layout_file = fullfile(script_dir, 'biosemi72.sfp');
        EEG = pop_chanedit(EEG, 'load', {layout_file, 'filetype', 'sfp'});

        % 2. Rereference
        EEG = pop_reref(EEG, []);

        % 3. High-pass filter
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.1);

        if run_ica_flag
            % 4. ICA (Apply 1Hz Rule and remove extreme samples)
            % Create copy for ICA and apply 1Hz high-pass
            EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
            
            % In EEGLAB, we must physically slice out bad continuous segments before ICA
            % (matches EegFun.jl's sample_selection logic - calculated on 0.1Hz data)
            bad_mask = any(abs(EEG.data) > 250, 1);
            bad_diff = diff([0, bad_mask, 0]);
            start_idx = find(bad_diff == 1);
            end_idx = find(bad_diff == -1) - 1;
            
            if ~isempty(start_idx)
                EEG_ica = pop_select(EEG_ica, 'nopoint', [start_idx(:), end_idx(:)]);
            end
            
            % Run extended ICA on the cleaned 1Hz filtered data
            % (Using 1e-5 to match the stringency of Julia's sum-of-squares 1e-6 stopping metric)
            EEG_ica = pop_runica(EEG_ica, 'extended', 1, 'interupt', 'off', 'pca', EEG_ica.nbchan - 1, 'stop', 1e-5);
            
            % Copy ICA weights back to the 0.1Hz filtered dataset
            EEG.icaweights = EEG_ica.icaweights;
            EEG.icasphere = EEG_ica.icasphere;
            EEG.icawinv = EEG_ica.icawinv;
            EEG.icachansind = EEG_ica.icachansind;
            
            % Sanity check: Plot the first 30 ICA components for the first subject
            if i == 1
                pop_topoplot(EEG_ica, 0, 1:30, 'EEGLAB ICA Components', [5 6]);
                f_ica = gcf;
                exportgraphics(f_ica, fullfile(data_dir, 'benchmark_matlab_ica.pdf'), 'ContentType', 'vector');
                close(f_ica);
            end
            
            % 5. Artifact Removal
            % Remove component 1
            EEG = pop_subcomp(EEG, 1, 0);
        end

        % 6. Sequence-based Trigger Processing
        % Convert numeric event types to strings so we can rename them to text labels safely
        for e = 1:length(EEG.event)
            if isnumeric(EEG.event(e).type)
                EEG.event(e).type = num2str(EEG.event(e).type);
            end
        end

        % We look for target event '5' preceded by valid (1, 3) or invalid (2, 4) cues.
        for e = 2:length(EEG.event)
            curr_type = EEG.event(e).type;
            prev_type = EEG.event(e-1).type;

            if strcmp(curr_type, '5')
                if strcmp(prev_type, '1') || strcmp(prev_type, '3')
                    EEG.event(e).type = 'valid_target';
                elseif strcmp(prev_type, '2') || strcmp(prev_type, '4')
                    EEG.event(e).type = 'invalid_target';
                end
            end
        end

        % 7. Epoch and Baseline
        % Apply 30Hz lowpass filter before epoching to avoid edge artifacts (matches EegFun.jl smoothing)
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 30.0);
        events = {'valid_target', 'invalid_target'};
        EEG = pop_epoch(EEG, events, [-0.5 1.0], 'newname', 'epochs', 'epochinfo', 'yes');
        EEG = pop_rmbase(EEG, [-200 0]);

        % Separate by condition to match EegFun's multi-condition extraction
        EEG_valid = pop_selectevent(EEG, 'type', 'valid_target');
        EEG_invalid = pop_selectevent(EEG, 'type', 'invalid_target');

        % Averaging into ERPs is usually handled by ERPLAB plugin or custom
        % code. We compute simple average matrices for benchmarking.
        all_erps_valid(:, :, i) = mean(EEG_valid.data, 3);
        all_erps_invalid(:, :, i) = mean(EEG_invalid.data, 3);
    end

    % 8. Grand Average
    grand_avg_valid = mean(all_erps_valid, 3);
    grand_avg_invalid = mean(all_erps_invalid, 3);

    % 10. Plot result (PO7 and PO8 average)
    po7_idx = find(strcmp({EEG.chanlocs.labels}, 'PO7'));
    po8_idx = find(strcmp({EEG.chanlocs.labels}, 'PO8'));

    avg_channel_erp_valid = mean(grand_avg_valid([po7_idx, po8_idx], :), 1);
    avg_channel_erp_invalid = mean(grand_avg_invalid([po7_idx, po8_idx], :), 1);

    times_sec = EEG.times(:) / 1000;
    export_table = table(times_sec, avg_channel_erp_valid(:), avg_channel_erp_invalid(:), 'VariableNames', {'time', 'valid', 'invalid'});
    writetable(export_table, fullfile(data_dir, 'benchmark_matlab_data.csv'));

    f = figure('Visible', 'off');
    plot(times_sec, avg_channel_erp_valid, 'LineWidth', 2, 'Color', 'b');
    hold on;
    plot(times_sec, avg_channel_erp_invalid, 'LineWidth', 2, 'Color', 'r');
    xline(0, 'k-', 'LineWidth', 1);
    yline(0, 'k-', 'LineWidth', 1);
    
    xlim([-0.2 1.0]);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title('Grand Average (PO7/PO8)');
    legend('Valid', 'Invalid', 'Location', 'northeast');
    grid off;

    exportgraphics(f, fullfile(data_dir, 'benchmark_matlab_erp.pdf'), 'ContentType', 'vector');
    
    t_elapsed = toc(t_start);
    
    % Save execution time
    fileID = fopen(fullfile(data_dir, 'benchmark_matlab_time.txt'), 'w');
    if fileID ~= -1
        fprintf(fileID, '%f', t_elapsed);
        fclose(fileID);
    end
    
    fprintf('EEGLAB execution time: %.2f seconds\n', t_elapsed);
end
