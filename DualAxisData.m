load('grad.mat')
n_chanels = length(grad.label);

%T = readtable('Sensor_positions_0mm.csv');

%% t

for chanel = n_chanels+1:n_chanels*2
    %labels
    label = grad.label(chanel-n_chanels);
    label = replace(label, 'bz', 'by');
    grad.label(chanel) = label;

    %chantype
    grad.chantype(chanel) = {'megmag'};

    %chanunit
    grad.chanunit(chanel) = {'T'};

    %chanpos/coilpos
    grad.chanpos(chanel,:) = grad.chanpos(chanel-n_chanels,:);
    grad.coilpos(chanel,:) = grad.coilpos(chanel-n_chanels,:);

    %chanori/coilori
    ori = T{(chanel-n_chanels+1), 10:12};
    ori_dot = strrep(ori, ',', '.');
    ori_double = str2double(ori_dot);
    grad.chanori(chanel,:) = ori_double;
end

%% Left
for i = 65:124
    ori = T{i, 10:12};
    ori_dot = strrep(ori, ',', '.');
    ori_double = str2double(ori_dot);
    row = i - 65 + n_chanels + 1;
    grad.chanori(row,:) = ori_double;
    grad.coilori(row,:) = ori_double;
end

%% Right

for i = 2:64
    ori = T{i, 10:12};
    ori_dot = strrep(ori, ',', '.');
    ori_double = str2double(ori_dot);
    row = 184 - 2 + i;
    grad.chanori(row,:) = ori_double;
    grad.coilori(row,:) = ori_double;
end

%% tra

grad.tra = eye(246);

%% Save
save('grad_dual_axis.mat', 'grad')

%% Save T
writetable(T, 'sensor_positions_123.csv');