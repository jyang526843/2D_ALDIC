function [poisson] = PlotPoisson0(DICpara,ResultStrain,sizeOfImg)
%PLOTPOISSON0: to compute and plot DIC solved Poisson's ratio 
%   [poisson]     = Plotstress0(DICpara,ResultStrain,sizeOfImg)
%
%   INPUT: DICpara          DIC para in the ALDIC code
%          ResultStrain     ALDIC computed strain field result
%          SizeOfImg        Size of the DIC raw image
%           
%   OUTPUT: poisson      Poisson's ratio
%
%   Plots:       
%       1) Poisson's ratio field
%       2) Poisson's ratio histogram
%
%   TODO: users could change caxis range based on their own choices.
%
% ----------------------------------------------
% Author: MFO  
% Contact and support: marcelo.falcao@usp.br
% Last time updated: 2023.11
% Reference: Plotstress0.m by Jin Yang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
um2px = DICpara.um2px; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)
   

%% Load computed strain fields
x2 = ResultStrain.strainxCoord; 
y2 = ResultStrain.strainyCoord;
%[M,N] = size(x2);
dudx = ResultStrain.dudx;
dvdy = ResultStrain.dvdy;
%% Compute Poisson's ratios
if dudx < 0 & dvdy > 0 
    poisson = -dudx./dvdy;
elseif dvdy < 0 & dudx > 0
    poisson = -dvdy./dudx;
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
[M,N] = size(x2);

%% Load displacement components to deform the reference image
% disp_u = ResultStrain.dispu; 
% disp_v = ResultStrain.dispv;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Poisson's ratio field ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  
surf(x2, um2px*(sizeOfImg(2)+1)-y2,poisson,'EdgeColor','none','LineStyle','none')
set(gca,'fontSize',18); view(2); box on; set(gca,'ydir','normal');
title('Poisson''s ratio','FontWeight','Normal','Interpreter','latex');
axis tight; axis equal; colorbar; colormap jet; set(gcf,'color','w');
if x2(M,N)<um2px*200, set(gca,'XTick',[]); end
if y2(M,N)<um2px*200, set(gca,'YTick',[]); end
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end

a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.0,0.01]);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

end
 
