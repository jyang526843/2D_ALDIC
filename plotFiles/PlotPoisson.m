function [poisson] = PlotPoisson(DICpara,ResultStrain,sizeOfImg,CurrentImg)
%PLOTPOISSON: to compute and plot DIC solved Poisson's ratio field on the original DIC images
%   [poisson]     = Plotstress(DICpara,ResultStrain,sizeOfImg,CurrentImg)
%
%   INPUT: DICpara          DIC para in the ALDIC code
%          ResultStrain     DIC computed strain field result
%          SizeOfImg        Size of the DIC raw image
%          CurrentImg       File name of current deformed image
%
%   OUTPUT: poisson         Poisson's ratio
%
%   Plots:       
%       1) Poisson's ratio field
%
%   TODO: users could change caxis range based on their own choices.
%
% ----------------------------------------------
% Author: MFO  
% Contact and support: 
% Last time updated: 2023.11
% Reference: Plotstress.m by Jin Yang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
warning off; load('colormap_RdYlBu.mat','cMap');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
um2px = DICpara.um2px; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)
   

%% Load computed strain fields
x2 = ResultStrain.strainxCoord; 
y2 = ResultStrain.strainyCoord;
dudx = ResultStrain.dudx; 
dvdy = ResultStrain.dvdy;

%% Compute Poisson's ratios
if dudx < 0 & dvdy > 0 | dvdy < 0 & dudx > 0
    if abs(dudx) < abs(dvdy)
        poisson = -dudx./dvdy;
    else 
        poisson = -dvdy./dudx;
    end
else
    fprintf ('ERROR! Deformation mismatch at some points.');
    return;
end
% trim borders, 5% each side
[hlen, vlen] = size(poisson);
trimh = int16(hlen*0.05);
if trimh < 1 
      trimh = 1;
end
trimv = int16(vlen*0.05);
if trimv < 1
     trimv = 1;
end
poisson = poisson(trimh+1:hlen-trimh,trimv+1:vlen-trimv);
x2 = x2(trimh+1:hlen-trimh,trimv+1:vlen-trimv);
y2 = y2(trimh+1:hlen-trimh,trimv+1:vlen-trimv);

%% Load displacement components to deform the reference image
disp_u = ResultStrain.dispu; 
disp_v = ResultStrain.dispv;
disp_u = disp_u(trimh+1:hlen-trimh,trimv+1:vlen-trimv);
disp_v = disp_v(trimh+1:hlen-trimh,trimv+1:vlen-trimv);
         
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Poisson plot field ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

title('Poisson''s ratio','FontWeight','Normal','Interpreter','latex');
axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=surf( (x2+Image2PlotResults*disp_u)/um2px, ...
    sizeOfImg(2)+1-(y2-Image2PlotResults*disp_v)/um2px, poisson,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; caxis auto; % set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight; colormap(cMap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.025,0.025]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) Poisson's ratio histogram ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
histogram(poisson,'FaceColor','white');  
xlabel('Poisson''s ratio','FontSize',12);
ylabel('Count','FontSize',12);
meanlabel = sprintf('Mean Poisson''s ratio: %.4f', mean(poisson(:)));
errorlabel = sprintf('Error (95%% CI): %.5f', 2*std(poisson(:)));
title({meanlabel,errorlabel});

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) Strain eyy histogram ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we trim 10% of eyy matrix borders
figure;
trimdvdy = dvdy(2*(trimh+1):hlen-2*trimh,2*(trimv+1):vlen-2*trimv);
histogram(trimdvdy,'FaceColor','white');  
xlabel('e_y_y','FontSize',12);
ylabel('Count','FontSize',12);
meanlabel = sprintf('Mean Strain e_y_y: %.5f', mean(trimdvdy(:)));
errorlabel = sprintf('Error (95%% CI): %.5f', 2*std(trimdvdy(:)));
title({meanlabel,errorlabel});

end
 
