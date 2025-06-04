import mne

folder_path_save = "/Users/lucke/Exjobb/MNEPython/"
file_prefix_save = "mne-signal-raw-"
dash = "-"
file_ext = ".fif"

folder_path = "/Users/lucke/Exjobb/Matlab/"
file_prefix = "signal-raw-"

n_signals = 10
projectors = 10

for signal in range(n_signals):
    for projector in range(projectors):
        path = f"{folder_path}{file_prefix}{signal+1}{file_ext}"
        raw = mne.io.read_raw_fif(path)
        empty_room_raw = mne.io.read_raw_fif('/Users/lucke/Exjobb/Matlab/signal-er-raw.fif')

        empty_room_projs = mne.compute_proj_raw(empty_room_raw, n_mag=projector+1)
        raw.add_proj(empty_room_projs, remove_existing=True)
        raw.apply_proj()

        file_name = f"{folder_path_save}{file_prefix_save}{signal+1}{dash}{projector+1}{file_ext}"
        raw.save(file_name, fmt='double')