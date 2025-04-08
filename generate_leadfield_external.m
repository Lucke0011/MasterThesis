function [lf_external, n_external_sources] = generate_leadfield_external()
% Load data structures
%load('grad.mat')
load('grad_dual_axis.mat')

%headmodel
cfg = [];
cfg.method = 'infinite';
headmodel_external = ft_prepare_headmodel(cfg);

%sourcemodel
ecg_pos = [0 0 -50];
env_1m_pos = [0 86.6 -50];
env_2m_pos = [-200 0 0];
env_3m_pos = [200 200 -100];
env_10m_pos = [0 1000 0];
dental_pos = [0 0.05 0.1];
pos_external = [ecg_pos; env_1m_pos; env_2m_pos; env_3m_pos; env_10m_pos; dental_pos];

cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = pos_external;
sourcemodel_external = ft_prepare_sourcemodel(cfg);

%leadfield
cfg = [];
cfg.grid           = sourcemodel_external;
cfg.headmodel      = headmodel_external;
cfg.grad           = grad;

leadfield_external_prepared = ft_prepare_leadfield(cfg); % compute leadfield
leadfield_external_prepared.leadfield = leadfield_external_prepared.leadfield';
n_external_sources = length(leadfield_external_prepared.leadfield);
n_sensors = length(leadfield_external_prepared.label);

% Current directions
ecg_vnn     = [1 0 0]; % Outwards
% Environment magnetic field direction perpendicular to vector direction
env_1m_vnn  = [0 1 0]; 
env_2m_vnn  = [1 0 0];
env_3m_vnn  = [1 1 0]/norm([1 1 0]);
env_10m_vnn = [0 1 0];
dental_vnn  = [0 1 0];
vnn_external = [ecg_vnn; env_1m_vnn; env_2m_vnn; env_3m_vnn; env_10m_vnn; dental_vnn];

lf_external = zeros(n_sensors, n_external_sources);
for i = 1:n_external_sources
    lf_external(:,i) = leadfield_external_prepared.leadfield{i} * vnn_external(i,:)';
end

lf_external = lf_external * 1e4;

end

