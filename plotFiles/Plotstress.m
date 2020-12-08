function [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
    stress_principal_min_xyplane, stress_maxshear_xyplane, ...
    stress_maxshear_xyz3d, stress_vonMises] = Plotstress(DICpara,ResultStrain,sizeOfImg,CurrentImg)
%PLOTSTRESS: to compute and plot DIC solved stress fields on the original DIC images
%   [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
%    stress_principal_min_xyplane, stress_maxshear_xyplane, ...
%    stress_maxshear_xyz3d, stress_vonMises]     = Plotstress(DICpara,ResultStrain,sizeOfImg,CurrentImg)
%
%   INPUT: DICpara          DIC para in the ALDIC code
%          ResultStrain     DIC computed strain field result
%          SizeOfImg        Size of the DIC raw image
%          CurrentImg       File name of current deformed image
%
%   OUTPUT: stress_sxx      Cauchy stress xx-compoent
%           stress_sxy      Cauchy stress xy-compoent s_xy = s_yx
%           stress_syy      Cauchy stress yy-compoent
%           stress_principal_max_xyplane    max principal stress on the xy-plane
%           stress_principal_min_xyplane    min principal stress on the xy-plane
%           stress_maxshear_xyplane         max shear stress on the xy-plane
%           stress_maxshear_xyz3d           max shear stress on the xyz-three dimensional space
%           stress_vonMises                 von Mises stress
%
%   Plots:       
%       1) stress sxx
%       2) stress sxy
%       3) stress syy
%       4) max principal stress on the xy-plane 
%       5) min principal stress on the xy-plane
%       6) max shear stress on the xy-plane
%       7) max shear stress on the xyz-three dimensional space
%       8) equivalent von Mises stress
%
% Author: Jin Yang  (jyang526@wisc.edu)
% Last date modified: 2020.11.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');
OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)
   
%% Load computed strain fields
x2 = ResultStrain.strainxCoord; 
y2 = ResultStrain.strainyCoord;
dudx = ResultStrain.dudx; dvdx = ResultStrain.dvdx; 
dudy = ResultStrain.dudy; dvdy = ResultStrain.dvdy; 
strain_exx = dudx; 
strain_exy = 0.5*(dvdx + dudy); 
strain_eyy = dvdy;

%% Load displacement components to deform the reference image
disp_u = ResultStrain.dispu; 
disp_v = ResultStrain.dispv;

%% Compute stress components

% ------ For each material model ------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DICpara.MaterialModel == 1 
    
    % Linear elasticity -- Plane stress
    % stress_xyz = [sxx  sxy  0;      strain_xyz = [exx  exy  0;
    %               sxy  syy  0;                    exy  eyy  0;
    %               0    0    0];                   0    0    ezz];
    E = DICpara.MaterialModelPara.YoungsModulus;
    nu = DICpara.MaterialModelPara.PoissonsRatio;
    stress_sxx = E/(1-nu^2)*(strain_exx + nu*strain_eyy);
    stress_syy = E/(1-nu^2)*(nu*strain_exx + strain_eyy);
    stress_sxy = E/(1+nu)*strain_exy;
     
    % Principal stress
    stress_maxshear_xyplane = sqrt((0.5*(stress_sxx-stress_syy)).^2 + stress_sxy.^2);
    stress_principal_max_xyplane = 0.5*(stress_sxx+stress_syy) + stress_maxshear_xyplane;
    stress_principal_min_xyplane = 0.5*(stress_sxx+stress_syy) - stress_maxshear_xyplane;
    stress_maxshear_xyz3d = max(stress_maxshear_xyplane, 0.5*abs(stress_principal_max_xyplane));
    
    % von Mises stress
    stress_vonMises = sqrt(0.5*( (stress_principal_max_xyplane-stress_principal_min_xyplane).^2 + ...
        (stress_principal_max_xyplane).^2 + (stress_principal_min_xyplane).^2 ));
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif DICpara.MaterialModel == 2 
    
    % Linear elasticity -- Plane strain
    % stress_xyz = [sxx  sxy  0;      strain_xyz = [exx  exy  0;
    %               sxy  syy  0;                    exy  eyy  0;
    %               0    0    szz];                 0    0    0]; 
    E = DICpara.MaterialModelPara.YoungsModulus;
    nu = DICpara.MaterialModelPara.PoissonsRatio;
    stress_sxx = E*(1-nu)/(1+nu)/(1-2*nu)*(strain_exx + nu/(1-nu)*strain_eyy);
    stress_syy = E*(1-nu)/(1+nu)/(1-2*nu)*(strain_eyy + nu/(1-nu)*strain_exx);
    stress_sxy = E/(1+nu)*strain_exy;
    stress_szz = nu*(stress_sxx + stress_syy);

    % Principal stress
    stress_maxshear_xyplane = sqrt((0.5*(stress_sxx-stress_syy)).^2 + stress_sxy.^2);
    stress_principal_max_xyplane = 0.5*(stress_sxx+stress_syy) + stress_maxshear_xyplane;
    stress_principal_min_xyplane = 0.5*(stress_sxx+stress_syy) - stress_maxshear_xyplane;
    
    stress_maxshear_xyz3d = reshape(  max( [ stress_maxshear_xyplane(:), 0.5*abs(stress_principal_max_xyplane(:)-stress_szz(:)), ...
                                0.5*abs(stress_principal_min_xyplane(:)-stress_szz(:)) ], [], 2 ),  size(stress_maxshear_xyplane) ) ;
    
     % von Mises stress
     stress_vonMises = sqrt(0.5*( (stress_principal_max_xyplane-stress_principal_min_xyplane).^2 + ...
        (stress_principal_max_xyplane-stress_szz).^2 + (stress_principal_min_xyplane-stress_szz).^2 ));
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif DICpara.MaterialModel == 3 % Neo-Hookean or other material models
    
    disp('User needs to modify the code by yourself.'); pause; % TODO: Please modify by yourself
    
end

 

         
%% ====== 1) Stress sxx ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+Image2PlotResults*disp_u,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),stress_sxx,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x2(end) < 200,set(gca,'XTick',[]); end
if y2(end) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('Stress $s_{xx}$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');



%% ====== 2) Strain sxy ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+Image2PlotResults*disp_u,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),stress_sxy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x2(end) < 200,set(gca,'XTick',[]); end
if y2(end) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('Stress $s_{xy}$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');



%% ====== 3) Strain syy ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+Image2PlotResults*disp_u,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),stress_syy,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.015,0.015]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x2(end) < 200,set(gca,'XTick',[]); end
if y2(end) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('Stress $s_{yy}$','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');

 

%% ====== 4) Strain stress_principal_max_xyplane ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+Image2PlotResults*disp_u,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),stress_principal_max_xyplane,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x2(end) < 200,set(gca,'XTick',[]); end
if y2(end) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$xy$-plane principal stress $s_{\max}$','FontWeight','Normal','Interpreter','latex');  set(gcf,'color','w');



%% ====== 5) Strain stress_principal_min_xyplane ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+Image2PlotResults*disp_u,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),stress_principal_min_xyplane,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x2(end) < 200,set(gca,'XTick',[]); end
if y2(end) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$xy$-plane principal stress $s_{\min}$','FontWeight','Normal','Interpreter','latex');  set(gcf,'color','w');



%% ====== 6) Strain stress_maxshear_xyplane ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+Image2PlotResults*disp_u,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),stress_maxshear_xyplane,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x2(end) < 200,set(gca,'XTick',[]); end
if y2(end) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$xy$-plane max shear stress','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');


%% ====== 7) Strain stress_maxshear_xyz3d ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+Image2PlotResults*disp_u,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),stress_maxshear_xyz3d,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x2(end) < 200,set(gca,'XTick',[]); end
if y2(end) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('$xyz$-3D max shear stress','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');

 
%% ====== 8) von Mises stress ======
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf(x2+Image2PlotResults*disp_u,sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v),stress_vonMises,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on;  caxis auto;  
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap jet; colormap(cMap);
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); % caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  %%Link them together
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = []; %%Hide the top axes
colormap(ax1,'gray'); % %%Give each one its own colormap

if x2(end) < 200,set(gca,'XTick',[]); end
if y2(end) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

xlabel( '$x$ (pixels)','Interpreter','latex');  ylabel('$y$ (pixels)','Interpreter','latex');
title('von Mises equivalent stress','FontWeight','Normal','Interpreter','latex'); set(gcf,'color','w');
 




end
 
 