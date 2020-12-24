function Plotdisp_show(U,coordinatesFEM,elementsFEM,varargin)
%FUNCTION Plotdisp_show(U,coordinatesFEM,elementsFEM)
% To plot DIC solved displacement components
% ----------------------------------------------
%
%   INPUT: U                 Displacement vector: 
%                            U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          coordinatesFEM    FE mesh coordinates
%          elementsFEM       FE mesh elements
%
%   OUTPUT: Plots of x-displacement field and y-displacement field.
%
%   TODO: users could change caxis range based on their own choices.
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 11/2020.
% ==============================================


%% Initialization
warning off;  U = full(U);

%% ====== 1) dispx u ======
figure; show([],elementsFEM(:,1:4),coordinatesFEM,U(1:2:end),'NoEdgeColor'); 

title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar; % view([90 -90])

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

colormap jet; box on;
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet);  % D Sample 
% caxis([-0.1,0.1]) % foam
% caxis([-0.004,0.004]); % Sample 12
% caxis([0,0.1]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ====== 2) dispy v ======
figure; show([],elementsFEM(:,1:4),coordinatesFEM,U(2:2:end),'NoEdgeColor'); 

title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

colormap jet; box on;
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet);  % D Sample 
% caxis([-0.1,0.1]) % foam
% caxis([-0.004,0.004]); % Sample 12
% caxis([0,0.1]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
