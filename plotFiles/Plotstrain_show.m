function Plotstrain_show(F,coordinatesFEM,elementsFEM,varargin)
%FUNCTION Plotstrain_show(F,coordinatesFEM,elementsFEM)
% To plot DIC solved strain components
% ----------------------------------------------
%
%   INPUT: F                 Deformation gradient tensor: 
%                            F = [F11_node1, F21_node1, F12_node1, F22_node1, ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
%          coordinatesFEM    FE mesh coordinates
%          elementsFEM       FE mesh elements
%
%   OUTPUT: Plots of exx, eyy, and exy strain fields.
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 11/2020.
% ==============================================

 
%% Initialization
warning off;  F = full(F);

%% ====== 1) strain exx ======
figure; show([],elementsFEM(:,1:4),coordinatesFEM,F(1:4:end),'NoEdgeColor'); 
title('Strain $e_{11}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal; 
colorbar;  box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;



%% ====== 2) strain exy ======
figure; show([],elementsFEM(:,1:4),coordinatesFEM,0.5*(F(2:4:end)+F(3:4:end)),'NoEdgeColor'); 
title('Strain $e_{12}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal;
colorbar;  box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;


%% ====== 3) strain eyy ======
figure; show([],elementsFEM(:,1:4),coordinatesFEM,F(4:4:end),'NoEdgeColor'); 
title('Strain $e_{22}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
view(2); axis tight; axis equal; 
colorbar;   box on;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
colormap jet;


