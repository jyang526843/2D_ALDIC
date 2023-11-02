% Save figures or output images for solved Poisson's ratio

%%
% Find img name
[~,imgname,imgext] = fileparts(file_name{1,ImgSeqNum});

%%
if DICpara.MethodToSaveFig == 1
    %% jpg
    figure(1); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis([-0.025,0.025]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_Poisson_ratio_field'],'-djpeg','-r300')
    
    figure(2);
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_Poisson_ratio_histogram'],'-djpeg','-r300')      
    
elseif DICpara.MethodToSaveFig == 2
    %% pdf
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_Poisson_ratio_field'];
    figure(1); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.025,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_Poisson_ratio_histogram'];
    figure(2);
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
   
else
    %% fig
    fprintf('Please modify codes manually in Section 8.');
    figure(1); colormap(coolwarm(128)); caxis([-0.05,0.1]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_Poisson_ratio_field']);
    figure(2); colormap(coolwarm(128)); caxis([-0.05,0.05]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_Poisson_ratio_histogram']);    
    
end