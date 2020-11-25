% ==============================================
% function Plotstrain_show
% ==============================================
 
function Plotstrain_show(F,coordinatesFEM,elementsFEM)

F = full(F);

figure; show([],elementsFEM(:,1:4),coordinatesFEM,F(1:4:end)); 
title('Strain $e_{11}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal; 
colorbar;  box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;


figure; show([],elementsFEM(:,1:4),coordinatesFEM,F(4:4:end)); 
title('Strain $e_{22}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal; 
colorbar;   box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;


figure; show([],elementsFEM(:,1:4),coordinatesFEM,0.5*(F(2:4:end)+F(3:4:end))); 
title('Strain $e_{12}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal;
colorbar;  box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;