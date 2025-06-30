%% Reset all
% clear all
% close all
% restoredefaultpath

%% Base paths
if contains(pwd,'/home/chrpfe')
    % Server:
    base_data_path = '/archive/23106_opmbci/';
    base_save_path = '/home/chrpfe/Documents/23106_opmbci';
    base_matlab_path = '/home/chrpfe/Documents/MATLAB/';
    project_scripts_path = '/home/chrpfe/Documents/MATLAB/23106_opmbci/opmbci_preprocessing';
else
    % Laptop:
    base_data_path = '/Volumes/KINGSTON/Recordings/Noise_recordings/MEG';
    base_save_path = '/Users/lucke/exjobb/opm_general-master';
    % base_matlab_path = '/Users/christophpfeiffer/Dropbox/Mac/Documents/MATLAB';
    project_scripts_path = '/Users/lucke/exjobb/opm_general-master';
end

load('CHOP251205_grad.mat')

%% Set up fieldtrip
% addpath(fullfile(base_matlab_path,'fieldtrip-20231220/')) % Fieldtrip path
% addpath(fullfile(base_matlab_path,'fieldtrip_private')) % Fieldtrip private functions
addpath(project_scripts_path)
ft_defaults

global ft_default
ft_default.showcallinfo = 'no';

%% Overwrite
overwrite = [];
overwrite.preproc = true;
overwrite.timelock = true;
overwrite.mri = false;
overwrite.coreg = false;
overwrite.sourcerec = false;

%% Params
params = [];
params.pre = 0.2; %sec
params.post = 0.8; %sec
params.pad = 0.2;
params.filter = [];
params.filter.hp_freq = 1;
params.filter.lp_freq = 70;
params.filter.bp_freq = [];
params.filter.notch = sort([50:50:150 60:60:120]);
params.n_comp = 40;
params.manual_ica = false;
params.ica_cor = 0.8; % cutoff for EOG/ECG coherence
params.ica_coh = 0.95; % cutoff for EOG/ECG coherence
params.save_ica = 1; % save plots and components
params.z_threshold = 20;
params.corr_threshold = 0.7; % correlation threshold for badchannel neighbors
params.opm_std_threshold = 2.5e-12;
params.hpi_freq = 33;
params.hpi_gof = 0.9;

% check trigger values
params.trigger_codes = [3 5 9];
params.trigger_labels = {'3', '5', '9'};

params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = '*bz';

%% Subjects + dates
sub = {'NatMEG_0953'};

ses = {'250203'};

%paradigms = {'RSEO'; 'RSEC'};
paradigms = {'RSEO'};

%%
i_sub = 1;
i_ses = 1;
%% Loop over subjects
params.sub = ['sub_' num2str(i_sub,'%02d')];
params.ses = ['ses_' num2str(i_ses,'%02d')];

%% Paths
raw_path = fullfile(base_data_path, sub{i_sub}, ses{i_ses});
save_path = fullfile(base_save_path, params.sub, params.ses);
if ~exist(save_path, 'dir')
   mkdir(save_path)
end
if ~exist(fullfile(save_path,'figs'), 'dir')
   mkdir(fullfile(save_path,'figs'))
end
for i_paradigm = 1:1
    %opm_files{i_paradigm} = fullfile(raw_path,'osmeg',[paradigms{i_paradigm} 'OPM_raw.fif']); % opm files 
    opm_files{i_paradigm} = '/Volumes/KINGSTON/Recordings/Noise_recordings/MEG/20241205_150002_sub-SBIRA27_file-VarITInoWirenoMove_raw.fif';
    aux_files{i_paradigm} = fullfile(raw_path,'meg',[paradigms{i_paradigm} 'EEG.fif']); % corresponding aux files containing EOG/ECG
end
hpi_path = fullfile(raw_path,'osmeg');
%mri_path = '/Volumes/dataarchvie/23106_opmbci/NatMEG_0953/mri/';
%mri_file = fullfile(mri_path,'mri','orig','001.mgz');

params.paradigm = paradigms{i_paradigm};

%% Read and preproc
if overwrite.preproc == true || ~exist(fullfile(save_path, [params.paradigm '_data_ica.mat']),'file')
    ft_hastoolbox('mne', 1);

    % Read data 
    disp(['Reading file: ' num2str(i_paradigm) '/' num2str(length(i_paradigm)) '...'])
    [data_epo, badchs_opm, opm_raw] = read_osMEG(opm_files{i_paradigm}, aux_files{i_paradigm}, save_path, params); % Read data
    
    save('data_epo.mat', 'data_epo');
    save('badchs_opm.mat', 'badchs_opm');

    %Reject visual here
    % cfg          = [];
    % cfg.method   = 'summary';
    % cfg.channel  = '*bz';
    % dummy        = ft_rejectvisual(cfg, data_epo);

    % Reject bad channels
    cfg = [];
    cfg.channel = setdiff(data_epo.label,badchs_opm);
    data_epo = ft_selectdata(cfg, data_epo);

    % Remove bad channels from grad
   

    % Reject jump trials
    % cfg = [];
    % cfg.channel = {'*bz'};
    % cfg.metric = 'maxzvalue';
    % cfg.preproc.medianfilter  = 'yes';
    % cfg.preproc.medianfiltord  = 9;
    % cfg.preproc.absdiff       = 'yes';
    % cfg.threshold = params.z_threshold;
    % [cfg,badtrl_jump] = ft_badsegment(cfg, data_epo);
    % data_epo = ft_rejectartifact(cfg,data_epo);
    
    % Reject noisy trials
    % cfg = [];
    % cfg.channel = {'*bz'};
    % cfg.metric = 'std';
    % cfg.threshold = params.opm_std_threshold;
    % [cfg,badtrl_std] = ft_badsegment(cfg, data_epo);
    % data_epo = ft_rejectartifact(cfg,data_epo);

    % Remove bad trials
    % [~,idx]=ismember(data_epo.sampleinfo,badtrl_jump,'rows');
    % badtrl_jump = find(idx);
    % [~,idx]=ismember(data_epo.sampleinfo,badtrl_std,'rows');
    % badtrl_std = find(idx);
    % save(fullfile(save_path, [params.paradigm '_badtrls']), ...
    %     'badtrl_jump', ...
    %     'badtrl_std', "-v7.3"); 
    
    % ICA
    % disp('Running ICA ...')
    % if sum(contains(data_ica.label,'EOG'))<1 || sum(contains(data_ica.label,'ECG'))<1 % No ExG data
    %     params.manual_ica = 1;
    %     params.save_ica = 1;
    % end
    % data_ica = ica_MEG(data_epo, save_path, params);
    % save(fullfile(save_path, [params.paradigm '_data_ica']), 'data_ica',"-v7.3"); disp('done');
end

% Remove bad trials
bad_trials = [4, 72, 138, 139, 166]';

cfg = [];
cfg.trials = setdiff(1:length(data_epo.trial), bad_trials);
data_epo = ft_selectdata(cfg, data_epo);

% Remove bad channels
bad_channels = {'ai57', 'ai58' 'ai59', 'ai60', 'hpiin271', 'hpiin272', 'hpiin273', 'hpiout271', 'hpiout272', 'hpiout273'};

cfg = [];
cfg.channel = setdiff(data_epo.label, bad_channels);
data_epo = ft_selectdata(cfg, data_epo);

% Demean
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-params.pre 0];
data_epo = ft_preprocessing(cfg, data_epo);

% Find NaN
% for i = 1:length(data_epo.trial)
%     nan_locs = isnan(data_epo.trial{i});
%     if any(nan_locs)
%         fprintf('nan in trial = %d\n', i);
%     end
% end

%% Mean data

mean_data = zeros(109, 1401);
for i = 1:length(data_epo.trial)
    mean_data = mean_data + data_epo.trial{i};
end

mean_data = mean_data / length(data_epo.trial);

%% HFC

orders = 3;
data_hfc = cell(1, orders); % signals x order
mean_data_hfc = cell(1, orders);
for i = 1:orders
    cfg = [];
    cfg.order = i;
    data_hfc{i} = ft_denoise_hfc(cfg, data_epo);

    temp = zeros(109, 1401);
    for j = 1:length(data_hfc{i}.trial)
        temp = temp + data_hfc{i}.trial{j};
    end
    
    mean_data_hfc{i} = temp / length(data_hfc{i}.trial);
end

%% Remove bad channels from grad
%grad_bad_channels = ['s91_bz', 's73_bz', 's32_bz', 'L114_bz'];
grad_bad_channels_i = [113, 112, 111, 13];

for i = 1:length(grad_bad_channels_i)
    grad.chanori(grad_bad_channels_i(i),:) = [];
    grad.chanpos(grad_bad_channels_i(i),:) = [];
    grad.chantype(grad_bad_channels_i(i),:) = [];
    grad.chanunit(grad_bad_channels_i(i),:) = [];
    grad.coilori(grad_bad_channels_i(i),:) = [];
    grad.coilpos(grad_bad_channels_i(i),:) = [];
    grad.label(grad_bad_channels_i(i),:) = [];
end
grad.tra = eye(109);

%%
data_epo.grad = grad;

%% Iterative SSS

data_sss = cell(1);
mean_data_sss = cell(1);
cfg = [];
cfg.sss.order_in = 8;
cfg.sss.order_out = 3;
data_sss{1} = ft_denoise_sss(cfg, data_epo);

temp = zeros(109, 1401);
for i = 1:length(data_sss{1}.trial)
    temp = temp + data_sss{1}.trial{i};
end

mean_data_sss{1} = temp / length(data_sss{1}.trial);

%% AMM

data_amm = cell(1);
mean_data_amm = cell(1);
cfg = [];
cfg.channel = 'all';
data_amm{1} = ft_denoise_amm(cfg, data_epo);

temp = zeros(109, 1401);
for i = 1:length(data_amm{1}.trial)
    temp = temp + data_amm{1}.trial{i};
end

mean_data_amm{1} = temp / length(data_amm{1}.trial);


%% Plot PSD

Fs = 1000;

%%

% HFC
figure;
subplot(2, 2, 1)
[psd, f] = pwelch(mean_data', [], [], 0:0.2:500, Fs);
semilogy(f, psd);
title('Signal Before HFC');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
xlim([1 70]);

subplot(2, 2, 2)
[psd, f] = pwelch(mean_data_hfc{1}', [], [], 0:0.2:500, Fs);
semilogy(f, psd);
title('Signal After HFC (1st Order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
xlim([1 70]);

subplot(2, 2, 3)
[psd, f] = pwelch(mean_data_hfc{2}', [], [], 0:0.2:500, Fs);
semilogy(f, psd);
title('Signal After HFC (2nd Order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
xlim([1 70]);

subplot(2, 2, 4)
[psd, f] = pwelch(mean_data_hfc{3}', [], [], 0:0.2:500, Fs);
semilogy(f, psd);
title('Signal After HFC (3rd Order)');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
xlim([1 70]);

%% Iterative SSS
figure;
subplot(2, 1, 1)
[psd, f] = pwelch(mean_data', [], [], 0:0.2:500, Fs);
semilogy(f, psd);
title('Signal Before Iterative SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
xlim([1 70]);

subplot(2, 1, 2)
[psd, f] = pwelch(mean_data_sss{1}', [], [], 0:0.2:500, Fs);
semilogy(f, psd);
title('Signal After Iterative SSS');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
xlim([1 70]);

%% AMM
figure
subplot(2, 1, 1)
[psd, f] = pwelch(mean_data', [], [], 0:0.2:500, Fs);
semilogy(f, psd);
title('Signal Before AMM');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
xlim([1 70]);

subplot(2, 1, 2)
[psd, f] = pwelch(mean_data_amm{1}', [], [], 0:0.2:500, Fs);
semilogy(f, psd);
title('Signal After AMM');
xlabel('Frequency (Hz)');
ylabel('PSD (T^2/Hz)');
grid on;
xlim([1 70]);

%% Butterfly plot

% label = 'hfc';
params.modality = 'opm';
params.layout = 'fieldlinebeta2bz_helmet.mat';
params.chs = '*bz';
params.amp_scaler = 1e15;
params.amp_label = 'B [fT]';
timelock(data_epo, params);
%end

%% clear and close all, then exit to free memory
close all
clear all
exit

%% Functions
function files = findOpmFiles(directory, pattern)
    % This function finds all files in the specified directory that match
    % the pattern "TrainingSet" followed by a number from 1 to 10 and then "_raw.fif".
    
    % Define the pattern
    pattern = '*_raw.fif';
    
    % Get a list of all files in the directory
    allFiles = dir(directory);
    
    % Initialize an empty cell array to store matching files
    files = {};
    
    % Loop through all files and check if they match the pattern
    for i = 1:length(allFiles)
        if ~allFiles(i).isdir
            if ~isempty(regexp(allFiles(i).name, pattern, 'once'))
                files{end+1} = fullfile(directory, allFiles(i).name); %#ok<AGROW>
            end
        end
    end
end