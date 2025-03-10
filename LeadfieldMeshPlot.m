%% Init

load('headmodels.mat')
load('sourcemodelT.mat')
load('grad.mat')

sourcemodelT.vnn = normals(sourcemodelT.pos,sourcemodelT.tri,'vertex');

%% generate the lead field

cfg = [];
cfg.grid = sourcemodelT;
cfg.grad = grad;
cfg.senstype         = 'meg'; 
cfg.headmodel = headmodels.headmodel_meg;
leadfield_prepared = ft_prepare_leadfield(cfg);
 
% convert leadfield from 3 orientations per source location to 1
% orientation (vnn) per source location and reshape to a matrix of size
% n_sensors x n_sources
LF = zeros(size(leadfield_prepared.leadfield{1},1),length(leadfield_prepared.leadfield));
for i = 1:length(leadfield_prepared.leadfield)
    LF(:,i) = leadfield_prepared.leadfield{i}*sourcemodelT.vnn(i,:)';
end

%%
% 1304
i_source = 9304; % picking a source

% plot the location and orientation of the source in the brain surface
figure
ft_plot_mesh(sourcemodelT, 'maskstyle', 'opacity', 'facecolor', 'black', 'edgecolor', 'red', 'unit','cm');
hold on
quiver3(sourcemodelT.pos(i_source,1),sourcemodelT.pos(i_source,2),sourcemodelT.pos(i_source,3),1e1*sourcemodelT.vnn(i_source,1),1e1*sourcemodelT.vnn(i_source,2),1e1*sourcemodelT.vnn(i_source,3),'LineWidth',4)
scatter3(sourcemodelT.pos(i_source,1),sourcemodelT.pos(i_source,2),sourcemodelT.pos(i_source,3), 5000,'b.')
hold off

%%
grad_mesh = plot_leadfield(grad.chanpos, LF(:,i_source));


% Function to make a mesh from the grad positions and plot with the edge
% color defined by the leadfield for one source (C = n_sensors x 1)
function mesh = plot_leadfield(points,C)
    % Check if the input is a Nx3 array
    if size(points, 2) ~= 3
        error('Input must be an Nx3 array of 3D points.');
    end

    % Calculate the center of the points
    center = mean(points, 1);

    % Translate points to center them around the origin
    centeredPoints = points - center;

    r = vecnorm(centeredPoints, 2, 2);
    theta = acos(centeredPoints(:, 3) ./ r);  % polar angle
    phi = atan2(centeredPoints(:, 1), centeredPoints(:, 2));  % azimuth angle with zero on y-axis
    x_proj = theta .* cos(phi);
    y_proj = theta .* sin(phi);

    % Perform Delaunay triangulation on the 2D polar coordinates
    tri = delaunay(x_proj, y_proj);

    % Create the mesh struct
    mesh.pos = points; % Original 3D points
    mesh.tri = tri;    % Triangulation indices

    % Display the results
    figure;
    trimesh(tri, points(:,1), points(:,2), points(:,3), C);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Triangulated Mesh');
    axis equal;
end

 