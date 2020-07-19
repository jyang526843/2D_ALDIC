function Plotdisp(U,x,y,sizeOfImg,CurrentImg,OrigDICImgTransparency)
  
load('./plotFiles/colormap_RdYlBu.mat','cMap');
%
% UWorld,x0,y0,size(ImgNormalized{1}),file_name{1,ImgSeqNum+1},OrigDICImgTransparency
% U = UWorld; x=x0;y=y0; sizeOfImg = size(ImgNormalized{1}); CurrentImg = file_name{1,ImgSeqNum+1}; OrigDICImgTransparency=0.5;

warning off;

M = size(x,1); N = size(x,2);
u = U(1:2:end); v = U(2:2:end);
u = reshape(u,M,N); v = reshape(v,M,N);

% imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
if M < 9 || (x(2)-x(1)<4), x2 = x(:,1)'; else x2 = interp(x(:,1)',4); end
if N < 9 || (y(1,2)-y(1,1)<4) ,y2 = y(1,:); else y2 = interp(y(1,:),4); end

z_u = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u,M*N,1),x2,y2);
z_v = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v,M*N,1),x2,y2);

[x2,y2]=ndgrid(x2,y2);x2=x2'; y2=y2';

%figure; 
% contourf(x2,y2,z,10,'EdgeColor','none','LineStyle','none')
% set(gca,'fontSize',18); view(2);
% set(gca,'ydir','normal');
% title('Strain $e_{xx}$','FontWeight','Normal','Interpreter','latex');
% axis tight; axis equal; colorbar; 
% if x(M,N) < 200
%     set(gca,'XTick',[]);
% end
% if y(M,N) < 200
%     set(gca,'YTick',[]);
% end
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
% 
% colormap jet
% box on

% alpha(c,.5)
fig1=figure; %1
ax1=axes;  
try
    h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch
    h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end
axis equal; axis tight; box on;
set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');

hold on; ax2=axes; 

h2=surf(x2+z_u ,sizeOfImg(2)+1-(y2-z_v),z_u,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency); axis equal;  axis tight;
colormap(coolwarm(32)); colormap jet; colormap(cMap);
%%Link them together
linkaxes([ax1,ax2]);  
%%Hide the top axes
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray');% set(ax2,'Colormap',coolwarmmap);

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
 
ax1.Visible = 'on';
ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); 
cb2.TickLabelInterpreter = 'latex';
set(gcf,'color','w');
xlabel( '$x$ (pixels)','Interpreter','latex'); 
ylabel('$y$ (pixels)','Interpreter','latex');
title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');






fig1=figure; %1
ax1=axes; 
try
    h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch
    h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end


axis on; axis equal; axis tight; box on;
set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');

hold on; ax2=axes; h2=surf(x2+z_u,sizeOfImg(2)+1-(y2-z_v),z_v,'EdgeColor','none','LineStyle','none');
set(gca,'fontSize',18); view(2); box on; %set(gca,'ydir','normal');
alpha(h2,OrigDICImgTransparency);  axis equal;  axis tight;
colormap(coolwarm(32));  colormap jet; colormap(cMap);
%%Link them together
linkaxes([ax1,ax2]);  
%%Hide the top axes
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray'); % colormap(ax2,coolwarmmap); 

if x(M,N) < 200,set(gca,'XTick',[]); end
if y(M,N) < 200,set(gca,'YTick',[]); end
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
 
ax1.Visible = 'on';
ax1.TickLabelInterpreter = 'latex'; 
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); 
cb2.TickLabelInterpreter = 'latex';
set(gcf,'color','w');
xlabel( '$x$ (pixels)','Interpreter','latex'); 
ylabel('$y$ (pixels)','Interpreter','latex');
title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');

 

