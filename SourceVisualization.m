vnn_out = [0 1 0];

% plot the location and orientation of the source
figure
%ft_plot_mesh(sourcemodel_out, 'maskstyle', 'opacity', 'facecolor', 'black', 'edgecolor', 'red', 'unit','cm');
scatter3(sourcemodel_out.pos(1,1),sourcemodel_out.pos(1,2),sourcemodel_out.pos(1,3), 5000,'b.')
hold on
quiver3(sourcemodel_out.pos(1,1),sourcemodel_out.pos(1,2),sourcemodel_out.pos(1,3),vnn_out(1,1),vnn_out(1,2),vnn_out(1,3),'LineWidth',4)
axis equal
title('Biological noise source and direction');
xlabel('X');
ylabel('Y');
zlabel('Z');
hold off