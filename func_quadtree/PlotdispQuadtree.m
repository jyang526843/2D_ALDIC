%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To plot DIC solved displacement components on the original DIC images
%   1) dispx 
function PlotdispQuadtree(U,coordinatesFEM,elementsFEM,CurrentImg,DICpara)
%FUNCTION PlotdispQuadtree(U,coordinatesFEMWorld,elementsFEM,CurrentImg,DICpara)
% To plot DIC solved displacement components  
% ----------------------------------------------
%
%   INPUT: U                    Displacement vector: U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          coordinatesFEM       FE mesh coordinates
%          elementsFEM          FE mesh elements
%          CurrentImg           Current deformed image
%          DICpara              DIC paramters
%
%   OUTPUT: Plots of x-displacement field and y-displacement field.
%
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 11/2020.
% ==============================================


%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');
OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

disp_u = U(1:2:end); disp_v = U(2:2:end);
coordinatesFEMWorldDef = [coordinatesFEM(:,1)+Image2PlotResults*disp_u, coordinatesFEM(:,2)+Image2PlotResults*disp_v];


%% ====== 1) dispx u ======
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; % h2=surf(x2+Image2PlotResults*disp_u ,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),disp_u,'EdgeColor','none','LineStyle','none');
 h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef,disp_u,'NoEdgeColor');
%h2=showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorldDef,disp_u);
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(coolwarm(32)); colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
colormap(jet);  
% caxis([-35,35]); % caxis([-0.025,0.025]); 
% caxis([-1.3,-0.1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap

if max(coordinatesFEM(:,1)) < 200,set(gca,'XTick',[]); end
if max(coordinatesFEM(:,2)) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');


 
%% ====== 2) dispy v ======
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; % h2=surf(x2+Image2PlotResults*disp_u ,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),disp_v,'EdgeColor','none','LineStyle','none');
h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef,disp_v,'NoEdgeColor');
%h2=showQuadtree(elementsFEM(:,1:4),coordinatesFEMWorldDef,disp_v);
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight; colormap(coolwarm(32)); colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
colormap(jet); 
% caxis([5,12]);
% caxis([0,38]); % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap

if max(coordinatesFEM(:,1)) < 200,set(gca,'XTick',[]); end
if max(coordinatesFEM(:,2)) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');




