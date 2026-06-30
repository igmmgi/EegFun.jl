import mne
import glob
import time
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg') # Non-interactive plotting
import matplotlib.pyplot as plt
from mne.preprocessing import ICA

def python_mne_benchmark(raw_files, data_dir, run_ica):

    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Dictionary to hold the single-subject evokeds for each condition
    all_subject_evokeds = {
        'valid': [], 'invalid': []
    }
    event_dict = {'valid': 101, 'invalid': 102}
    
    results_dir = os.path.join(data_dir, "benchmarks")
    os.makedirs(results_dir, exist_ok=True)
    
    for i, file in enumerate(raw_files):

        # 1. Read data and layout file
        raw = mne.io.read_raw_bdf(file, preload=True)
        montage_path = os.path.join(script_dir, "biosemi72.sfp")
        montage = mne.channels.read_custom_montage(montage_path)
        raw.set_montage(montage, match_case=False)
        
        # 72 EEG Channel + 1 Status Channel (BDF data)
        ch_types = {ch: 'eeg' for ch in raw.ch_names if ch != 'Status'}
        raw.set_channel_types(ch_types)
        
        # 2. Initial preprocessing options: rereference, high-pass filter, and mark bad data sections
        raw.set_eeg_reference('average')
        # 2. High-pass filter
        raw.filter(l_freq=0.1, h_freq=None, phase='zero')

        data_raw = raw.get_data(picks='eeg')
        bad_mask = np.any(np.abs(data_raw) > 250e-6, axis=0)
        good_mask = ~bad_mask

        if run_ica:
            # 3. ICA (Apply 1Hz Rule and remove extreme samples)
            raw_ica = raw.copy().filter(l_freq=1.0, h_freq=None, phase='zero')
            data_ica_full = raw_ica.get_data()
            clean_raw_ica = mne.io.RawArray(data_ica_full[:, good_mask], raw_ica.info)
            
            n_eeg = len(mne.pick_types(raw.info, eeg=True, exclude=[]))
            ica = ICA(n_components=n_eeg - 1, method='infomax', fit_params=dict(extended=True, verbose=True))
            ica.fit(clean_raw_ica)
                 
            # Sanity check: Plot the first 30 ICA components for the first subject
            if i == 0:
                fig_ica = ica.plot_components(picks=range(30), show=False)
                fig_ica.savefig(os.path.join(results_dir, "python_ica.pdf"))
                plt.close('all')
                 
            # 4. Artifact Removal
            ica.exclude = [0]
            ica.apply(raw)
        
        # 5. Sequence-based Trigger Processing
        events = mne.find_events(raw)
        is_target = events[1:, 2] == 5
        prev_valid = np.isin(events[:-1, 2], [1, 3])
        prev_invalid = np.isin(events[:-1, 2], [2, 4])
        
        valid_events = events[1:][is_target & prev_valid].copy()
        valid_events[:, 2] = 101  # Recode as 'valid'
        invalid_events = events[1:][is_target & prev_invalid].copy()
        invalid_events[:, 2] = 102  # Recode as 'invalid'
        
        new_events = np.vstack((valid_events, invalid_events))
        
        # 6. Epoch and Baseline
        epochs = mne.Epochs(
            raw, new_events, event_dict, 
            tmin=-0.5, tmax=1.0, 
            baseline=(-0.2, 0.0), 
            preload=True
        )
        
        # Average trials and apply 30Hz lowpass filter (matches EegFun.jl)
        for cond in event_dict.keys():
            evoked = epochs[cond].average()
            # 7. Low-pass filter on the ERP (Moved to Evoked level to match Julia pipeline order)
            evoked.filter(l_freq=None, h_freq=30.0, phase='zero')
            all_subject_evokeds[cond].append(evoked)
                
    # 7. Grand Average
    grand_averages = {
        cond: mne.grand_average(evokeds) 
        for cond, evokeds in all_subject_evokeds.items()
    }
    
    # 8. Plot result (PO7 and PO8 average)
    fig, ax = plt.subplots()
    picks = mne.pick_channels(grand_averages['valid'].info['ch_names'], ['PO7', 'PO8'])
    
    valid_data = grand_averages['valid'].data[picks, :].mean(axis=0) * 1e6 # Convert to uV
    invalid_data = grand_averages['invalid'].data[picks, :].mean(axis=0) * 1e6
    times = grand_averages['valid'].times # Keep in seconds
    
    export_data = np.column_stack((times, valid_data, invalid_data))
    np.savetxt(os.path.join(results_dir, "python_data.csv"), export_data, delimiter=",", header="time,valid,invalid", comments="")

    
    ax.plot(times, valid_data, label='Valid', color='b', linewidth=2)
    ax.plot(times, invalid_data, label='Invalid', color='r', linewidth=2)
    ax.axvline(0, color='k', linewidth=1)
    ax.axhline(0, color='k', linewidth=1)
    ax.set_xlim(-0.2, 1.0)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude (µV)')
    ax.set_title('Grand Average (PO7/PO8)')
    ax.legend(loc='upper right')
    ax.grid(False)
    
    fig.savefig(os.path.join(results_dir, "python_erp.pdf"))
    plt.close(fig)
    
    return grand_averages

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if len(sys.argv) < 2:
        raise ValueError("Please provide the data directory as a command line argument (e.g., python mne_benchmark.py /path/to/data)")
    
    data_dir = sys.argv[1]
        
    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"Data directory not found: {data_dir}")
    
    test_files = sorted(glob.glob(os.path.join(data_dir, "*.bdf")))
    
    n_files_to_process = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    if n_files_to_process > 0:
        test_files = test_files[:n_files_to_process]
        
    run_ica_flag = sys.argv[3].lower() in ['true', '1', 't', 'y', 'yes'] if len(sys.argv) > 3 else True
    
    print("Benchmarking MNE-Python pipeline...")
    t0 = time.time()
    python_mne_benchmark(test_files, data_dir, run_ica_flag)
    t1 = time.time()
    
    elapsed = t1 - t0
    results_dir = os.path.join(data_dir, "benchmarks")
    with open(os.path.join(results_dir, "python_time.txt"), "w") as f:
        f.write(str(elapsed))
        
    print(f"MNE-Python execution time: {elapsed:.2f} seconds")
