function [lf_brain, n_brain_sources] = generate_leadfield_brain()
    % Load data structures
    load('sourcemodelT.mat')
    load('headmodels.mat')
    %load('grad.mat')
    load('grad_dual_axis.mat')
    
    % Calculate dipole directions
    sourcemodelT.vnn = normals(sourcemodelT.pos,sourcemodelT.tri,'vertex');
    
    % Leadfield
    cfg.grid           = sourcemodelT;
    cfg.headmodel      = headmodels.headmodel_meg;
    cfg.grad           = grad;
    leadfield_prepared = ft_prepare_leadfield(cfg);
    n_brain_sources = length(leadfield_prepared.leadfield);
    n_sensors = length(leadfield_prepared.label);
    
    lf_brain = zeros(n_sensors,n_brain_sources);
    for i = 1:n_brain_sources
        vnnT = sourcemodelT.vnn(i,:)';
        lf_brain(:,i) = leadfield_prepared.leadfield{i} * vnnT;
    end
    
    lf_brain = lf_brain * 1e4;

end

