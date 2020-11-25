% ---------------------------------------------
% Augmented Lagrangian DIC
% Author: Jin Yang, PhD. Caltech
% Contact: yangjin@caltech.edu
% 2015.04,06,07; 2016.03,04; 2020.11
% ---------------------------------------------

%% Section 1: Clear MATLAB environment & mex set up Spline interpolation  
close all; clear; clc; clearvars -global
fprintf('------------ Section 1 Start ------------ \n')
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
% % cd("./Splines_interp/lib_matlab"); CompileLib; cd("../../");  % % mex bi-cubic spline interpolations
% % addpath("./Splines_interp/lib_matlab"); % dbstop if error % % Old version codes.
mex -O ba_interp2.cpp; 
% Comment: If your MATLAB has just been crashed, and this line reports error but it works before the crash, 
% Change this line to: "try mex -O ba_interp2.cpp; catch; end"
addpath("./func",'./src','./plotFiles','./func_quadtree','./func_quadtree/refinement','./plotFiles/export_fig-d966721/'); 
% Also: addpath("./YOUR IMAGE FOLDER"); 
fprintf('------------ Section 1 Done ------------ \n \n')


%% Section 2: Load DIC parameters and set up DIC parameters 
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ====== 
[file_name,Img,DICpara] = ReadImageQuadtree; close all;
% %%%%%% Uncomment the line below to change the DIC computing region (ROI) manually %%%%%%
%gridxROIRange = [gridxROIRange1,gridxROIRange2]; gridyROIRange = [Val1, Val2];
%gridxROIRange = [224,918]; gridyROIRange = [787,1162];

% ====== Normalize images: fNormalized = (f-f_avg)/(f_std) ======
[ImgNormalized,DICpara.gridxyROIRange] = funNormalizeImg(Img,DICpara.gridxyROIRange); 

% ====== Initialize variable storage ======
ResultDisp = cell(length(ImgNormalized)-1,1);   ResultDefGrad = cell(length(ImgNormalized)-1,1);
ResultStrain = cell(length(ImgNormalized)-1,1); ResultStress = cell(length(ImgNormalized)-1,1);
ResultFEMesh = cell(ceil((length(ImgNormalized)-1)/DICpara.ImgSeqIncUnit),1); % For incremental DIC mode
fprintf('------------ Section 2 Done ------------ \n \n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fNormalized = ImgNormalized{1}; 
for ImgSeqNum = 2:length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update an initial guess of displacements
    % The key idea is to either to use a new FFT-based cross correlation peak fitting, 
    % or use the results from the last frame as the new initial guess for the next frame;
    % Particularly in incremental mode DIC, the reference image can also be updated, e.g.:
    % fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)};
    gNormalized = ImgNormalized{ImgSeqNum}; NewFFTSearchCheck = 0; DICpara.NewFFTSearch = 0;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1
        % ====== Integer Search ======
        [DICpara,x0temp,y0temp,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara);
        % ====== FEM mesh set up ======
        [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp;
        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0); %PlotuvInit; [x0temp,y0temp,u,v,cc]= IntegerSearchMg(fNormalized,gNormalized,file_name,DICpara);
        % ====== Deal with incremental mode ======
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1, fNormalizedNewIndex = fNormalizedNewIndex-1; end
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
       
    elseif mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0 % TO update ref image in incremental mode
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1,  fNormalizedNewIndex = fNormalizedNewIndex-1; end
        fNormalized = ImgNormalized{fNormalizedNewIndex}; % Update reference
        [DICpara,DICmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DICpara,DICmesh); % Update reference image if needed;
        U0 = zeros(2*size(DICmesh.coordinatesFEM,1),1); % PlotuvInit;
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
    else
        U0 = ResultDisp{ImgSeqNum-2}.U;
    end
    %Plotdisp_show(U0,[coordinatesFEM(:,1),size(fNormalized,2)+1-coordinatesFEM(:,2)],elementsFEM); % Plot initial values
    % ====== Spline interpolation images ======
    %[imgfNormalizedbc,imggNormalizedbc,imgSize,DfAxis] = funImgGradient(fNormalized,gNormalized);
    Df = funImgGradient(fNormalized,gNormalized); % % using finite difference;
    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    fprintf('------------ Section 3 Done ------------ \n \n')

    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section 3*: Generate an image mask and remove finite elements where there is a hole
    % This section is to deal with the hole geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1
        % ====== Generate an image mask in the reference image ======
        GenerateImageMask;
        % ====== Generate Quadtree mesh ======
        GenerateQuadtreeMesh;
    end
    
    %% Section 4: Subproblem 1 -or- Local ICGN Subset DIC
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve first local step in ALDIC: Subproblem 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu = 0; beta = 0; tol = 1e-2; ALSolveStep = 1; ALSub1Time = zeros(6,1); ALSub2Time = zeros(6,1); 
    ConvItPerEle = zeros(size(DICmesh.coordinatesFEM,1),6); ALSub1BadPtNum = zeros(6,1);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = ...
        LocalICGN(U0,DICmesh.coordinatesFEM,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp; toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    % Comment these lines below if you don't have local bad points
    % Comment START
    [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,USubpb1,FSubpb1,0);
    disp('--- Remove bad points done ---')
    % Comment END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ------ Plot ------
    USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); 
    FSubpb1World = FSubpb1; % FSubpb1World(2:2:end) = -FSubpb1(2:2:end);
    % close all; Plotdisp_show(USubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
    % Plotstrain_show(FSubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    fprintf('------------ Section 4 Done ------------ \n \n')

   
    %% Section 5: Subproblem 2 -or- solve the global compatible displacement field
    fprintf('------------ Section 5 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve global step in ALDIC: Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; LevelNo=1;
    % FSubpb1 = funSmoothStrainQuadtree(FSubpb1,DICmesh,DICpara);
    
	% ====== Define penalty parameter ======
    mu = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1; 
    betaList = [1e-3]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
    %Subpb2FDOrFEM = Using FE method for quadtree mesh
    M = size(DICmesh.x0,1); N = size(DICmesh.x0,2); GaussPtOrder = 2; alpha = 0;
    close all; % hbar = waitbar(0,'Please wait for Subproblem 2 global step!');
    % ====== Solver using finite element method ======
    for tempk = 1:length(betaList)
        display(['Try #',num2str(tempk),' beta = ',num2str(betaList(tempk))]);
        beta = betaList(tempk);
        GaussPtOrder=2; alpha=0; [USubpb2] = Subpb2Quadtree(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        % Plotdisp_show(USubpb2,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4));
        [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,GaussPtOrder,0);
        % Plotstrain_show(FSubpb1,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4));
          
        Err1(tempk) = norm(USubpb1-USubpb2,2);
        Err2(tempk) = norm(FSubpb1-FSubpb2,2);
        % waitbar(tempk/(length(betaList)+1));
    end
    Err1Norm = (Err1-mean(Err1))/std(Err1); figure, plot(Err1Norm);
    Err2Norm = (Err2-mean(Err2))/std(Err2); figure, plot(Err2Norm);
    ErrSum = Err1Norm+Err2Norm; figure,plot(ErrSum); title('Tune the best \beta in the subproblem 2'); [~,indexOfbeta] = min(ErrSum); 
    try % Tune the best beta by a quadratic polynomial fitting
        [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
        p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
    catch
        beta = betaList(indexOfbeta);
    end
    display(['Best beta = ',num2str(beta)]);
    % Using the optimal beta to solve again
    if abs(beta-betaList(end))>abs(eps)
        [USubpb2] = Subpb2Quadtree(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        USubpb2 = full(USubpb2); waitbar(1); close(hbar); ALSub2Time(ALSolveStep) = toc; toc
        % ------- Before computing strain, we smooth the displacement field -------
        % USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh,DICpara);
        % ------- Compute strain field --------
        [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,GaussPtOrder,0);
    end
    for tempk=0:3
    	FSubpb2(4*markCoordHoleEdge-tempk) = FSubpb1(4*markCoordHoleEdge-tempk);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------- Smooth strain field --------
    % FSubpb2 = funSmoothStrainQuadtree(FSubpb2,DICmesh,DICpara);

    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end); FSubpb2World = FSubpb2;  
    close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
    Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));

    % ======= Update dual variables =======
    udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
	save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Section 6: ADMM iterations
    fprintf('------------ Section 6 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to run ADMM iteration: Subproblem 1 & 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ==================== ADMM AL Loop ==========================
    ALSolveStep = 1; tol2 = 1e-4; UpdateY = 1e4; CrackOrNot = 0; CrackPath1 = [0,0]; CrackPath2 = [0,0]; CrackTip = [0,0]; 
    HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end

    while (ALSolveStep < 3)
        ALSolveStep = ALSolveStep + 1;  % Update using the last step
        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic;[USubpb1,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1(USubpb2,FSubpb2,udual,vdual,DICmesh.coordinatesFEM,...
        Df,fNormalized,gNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        FSubpb1 = FSubpb2; toc 
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------  Manually find some bad points from Local Subset ICGN step ------
        % disp('--- Start to manually remove bad points --- \n')
        % disp('    Comment codes here if you do not have bad local points \n')
        % % Comment START
        % [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,USubpb1,FSubpb1,0);
        % % Plotdisp_show(USubpb1,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4));
        % disp('--- Remove bad points done ---')
        % Comment END
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
        % USubpb1 = funSmoothDispQuadtree(USubpb1,DICmesh,DICpara);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ============== Subproblem 2 ==============
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        tic; [USubpb2] = Subpb2Quadtree(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
		ALSub2Time(ALSolveStep) = toc; toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ------- Before computing strain, we smooth the displacement field -------
        % USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh,DICpara);
        GaussPtOrder = 2; [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,GaussPtOrder,0);
        for tempk=0:3
            FSubpb2(4*markCoordHoleEdge-tempk) = FSubpb1(4*markCoordHoleEdge-tempk);
        end
         
		% ------- Smooth strain field --------
        % FSubpb2 = funSmoothStrainQuadtree(FSubpb2,DICmesh,DICpara);
        save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute norm of UpdateY
        USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
        USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
        USubpb1_Old = load(['Subpb1_step',num2str(ALSolveStep-1)],'USubpb1');
        USubpb1_New = load(['Subpb1_step',num2str(ALSolveStep)],'USubpb1');
        if (mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) ~= 0 && (ImgSeqNum>2)) || (ImgSeqNum < DICpara.ImgSeqIncUnit)
            UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(size(USubpb2_Old.USubpb2,1));
            try
                UpdateY2 = norm((USubpb1_Old.USubpb1 - USubpb1_New.USubpb1), 2)/sqrt(size(USubpb1_Old.USubpb1,1));
            catch
            end
        end
        try
            disp(['Update local step  = ',num2str(UpdateY2)]);
            disp(['Update global step = ',num2str(UpdateY)]);
        catch
        end
        fprintf('*********************************** \n \n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update dual variables------------------------------
        udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1; 
		 
        save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
        try
        if UpdateY < tol2 || UpdateY2 < tol2
            break
        end
        catch
        end

    end
    fprintf('------------ Section 6 Done ------------ \n \n')
 
    % Save data
    ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
    ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
    ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2); % tempFoamAL;
    
end    


%% ------ Plot ------
% ------- Smooth strain field --------
FSubpb2 = funSmoothStrainQuadtree(FSubpb2,DICmesh,DICpara);
USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
FSubpb2World = FSubpb2; % FSubpb1World(2:2:end) = -FSubpb1(2:2:end);
close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));

% ------ Save results ------
% Find img name and save all the results 
[~,imgname,imgext] = fileparts(file_name{1,end});
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ALSub1Time','ALSub2Time','ALSolveStep');

 

%% Section 7: Check convergence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to check convergence of ADMM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('------------ Section 7 Start ------------ \n')
% % ====== Check convergence ======
% fprintf('***** Check convergence ***** \n');
% ALSolveStep1 = min(6,ALSolveStep);
% disp('========');
% for ALSolveStep = 1:ALSolveStep1
%     USubpb2 = load(['Subpb2_step',num2str(ALSolveStep )],'USubpb2');
%     USubpb1 = load(['Subpb1_step',num2str(ALSolveStep )],'USubpb1');
%     UpdateY = norm((USubpb2.USubpb2 - USubpb1.USubpb1), 2)/sqrt(length(USubpb2.USubpb2));
%     disp(num2str(UpdateY));
% end
% disp('========');
% for ALSolveStep = 1:ALSolveStep1
%     FSubpb1 = load(['Subpb1_step',num2str(ALSolveStep )],'FSubpb1');
%     FSubpb2 = load(['Subpb2_step',num2str(ALSolveStep )],'FSubpb2');
%     UpdateF = norm((FSubpb1.FSubpb1 - FSubpb2.FSubpb2), 2)/sqrt(length(FSubpb1.FSubpb1));
%     disp(num2str(UpdateF));
% end
% disp('========');
% for ALSolveStep = 2:ALSolveStep1
%     USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
%     USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
%     UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(length(USubpb2.USubpb2));
%     disp(num2str(UpdateY));
% end
% disp('========');
% for ALSolveStep = 2:ALSolveStep1
%     uvdual_Old = load(['uvdual_step',num2str(ALSolveStep-1)],'udual');
%     uvdual_New = load(['uvdual_step',num2str(ALSolveStep)],'udual');
%     UpdateW = norm((uvdual_Old.udual - uvdual_New.udual), 2)/sqrt(length(uvdual_Old.udual));
%     disp(num2str(UpdateW));
% end
% disp('========');
% for ALSolveStep = 2:ALSolveStep1
%     uvdual_Old = load(['uvdual_step',num2str(ALSolveStep-1)],'vdual');
%     uvdual_New = load(['uvdual_step',num2str(ALSolveStep)],'vdual');
%     Updatev = norm((uvdual_Old.vdual - uvdual_New.vdual), 2)/sqrt(length(uvdual_Old.vdual));
%     disp(num2str(Updatev));
% end
% fprintf('------------ Section 7 Done ------------ \n \n')
% ------ Delete temp files ------
for tempi = 1:ALSolveStep
    file_name_Subpb1 = ['Subpb1_step',num2str(tempi),'.mat'];
    file_name_Subpb2 = ['Subpb2_step',num2str(tempi),'.mat'];
    file_name_dual = ['uvdual_step',num2str(tempi),'.mat'];
    delete(file_name_Subpb1); delete(file_name_Subpb2); delete(file_name_dual);
end

% ------ clear temp variables ------
clear a ALSub1BadPtNum ALSub1Timetemp atemp b btemp cc ConvItPerEletemp hbar Hbar 



%% Section 8: Compute strains
fprintf('------------ Section 8 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain fields and plot figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Smooth displacements ------
DICpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');                 
% ------ Choose strain computation method ------
DICpara.MethodToComputeStrain = 3; %funParaInput('StrainMethodOp'); 
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DICpara.StrainType = funParaInput('StrainType');
% ------ Choose image to plot (first only, second and next images) ------
if length(ImgNormalized)==2, DICpara.Image2PlotResults = funParaInput('Image2PlotResults');
else DICpara.Image2PlotResults = 1; % Plot over current, deformed image by default
end
% ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat');
% ------ Choose overlay image transparency ------
DICpara.OrigDICImgTransparency = 1; 
if DICpara.MethodToSaveFig == 1  
    DICpara.OrigDICImgTransparency = funParaInput('OrigDICImgTransparency');         
end

% ------ Start main part ------
for ImgSeqNum = 2:length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
    if DICpara.ImgSeqIncUnit > 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit);
    elseif DICpara.ImgSeqIncUnit == 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)-1;
    end
    FEMeshInd = FEMeshIndLast + 1;
    
    if FEMeshInd == 1
        USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
        coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
        elementsFEM = ResultFEMesh{1}.elementsFEM;
        if (ImgSeqNum-1 == 1) || (DICpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    else
        USubpb2 = ResultDisp{ImgSeqNum-1}.U;
        if mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0
            coordinatesFEM = ResultFEMesh{FEMeshInd}.coordinatesFEM;
            elementsFEM = ResultFEMesh{FEMeshInd}.elementsFEM;
            coordinatesFEMLast = ResultFEMesh{FEMeshIndLast}.coordinatesFEM;
            UFEMeshLast = ResultDisp{ImgSeqNum-2}.U + UFEMesh;
            xq = coordinatesFEM(:,1); yq = coordinatesFEM(:,2);
            UFEMesh = 0*USubpb2;
            UFEMesh(1:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(1:2:end),xq,yq,'v4');
            UFEMesh(2:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(2:2:end),xq,yq,'v4');
        end
        USubpb2 = USubpb2 + UFEMesh;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    coordinatesFEMWorld = [coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------ Plotting and Compute Strain-------
    if size(USubpb2,1) == 1
        ULocal = USubpb2_New.USubpb2; FLocal = FSubpb2.FSubpb2; 
    else
        ULocal = USubpb2; FLocal = FSubpb2;
    end
    % ULocal(1:2:end)= -ULocal(1:2:end); 
    % FLocal(1:4:end)= -1*FLocal(1:4:end); FLocal(3:4:end)= -1*FLocal(3:4:end); % because u is flipped sign
    UWorld = ULocal; UWorld(2:2:end) = -UWorld(2:2:end); FWorld = FLocal; % close all; Plotuv(UWorld,x0,y0World);
     
    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            ULocal = funSmoothDispQuadtree(ULocal,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
    
    % ----- Compute strain field ------
    ComputeStrainQuadtree; % Compute strain
    % % ------- Add filter and plot strain field -------
    %Plotstrain_Fij;

    % ------ Plot disp and strain ------
    close all;
    if DICpara.OrigDICImgTransparency == 2
         Plotdisp_show(UWorld,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
        [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = ...
                   Plotstrain0Quadtree(FStraintemp,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
            PlotdispQuadtree(UWorld,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),file_name{1,1},DICpara);
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
                strain_maxshear,strain_vonMises] = PlotstrainQuadtree(UWorld,FStraintemp, ...
                DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),file_name{1,1},DICpara);

        else % Plot over second or next deformed images
            PlotdispQuadtree(UWorld,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),file_name{1,ImgSeqNum},DICpara);
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
                strain_maxshear,strain_vonMises] = PlotstrainQuadtree(UWorld,FStraintemp, ...
                DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),file_name{1,ImgSeqNum},DICpara);
        
        end
    end

    % ----- Save strain results ------
    ResultStrain{ImgSeqNum-1} = struct('strainxCoord',DICmesh.coordinatesFEM(:,1),'strainyCoord',DICmesh.coordinatesFEM(:,2), ...
            'dispu',ULocal(1:2:end),'dispv',ULocal(2:2:end), ...
            'dudx',FStraintemp(1:4:end),'dvdx',FStraintemp(2:4:end),'dudy',FStraintemp(3:4:end),'dvdy',FStraintemp(4:4:end), ...
            'strain_exx',strain_exx,'strain_exy',strain_exy,'strain_eyy',strain_eyy, ...
            'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min, ...
            'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises);
     
    % ------ Save figures for tracked displacement and strain fields ------
    SaveFigFilesDispAndStrainQuadtree;
    
    
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 8 Done ------------ \n \n')

 
% ------ Save data again including solved strain fields ------
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh',...
                   'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain');



%% Section 9: Compute stress
fprintf('------------ Section 9 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute stress fields and plot figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Choose material model ------ 
DICpara.MaterialModel = funParaInput('MaterialModel');
% ------ Define parameters in material models ------
if (DICpara.MaterialModel == 1) || (DICpara.MaterialModel == 2) % Linear elasticity
    fprintf('Define Linear elasticity parameters \n')
    fprintf("Young's modulus (unit: Pa). "); prompt = 'Input here: '; 
    DICpara.MaterialModelPara.YoungsModulus = input(prompt); 
    fprintf("Poisson's ratio. "); prompt = 'Input here: '; 
    DICpara.MaterialModelPara.PoissonsRatio = input(prompt);
    fprintf('------------------------------------- \n');
end

% ------ Start main part ------
for ImgSeqNum = 2:length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]); close all;
     
    % ------ Plot disp and strain ------
    if DICpara.OrigDICImgTransparency == 2
        [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises]  =  Plotstress0Quadtree( ...
            DICpara,ResultStrain{ImgSeqNum-1},DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4)); 
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
            [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = PlotstressQuadtree( ...
                DICpara,ResultStrain{ImgSeqNum-1},DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),file_name{1,1});
             
        else % Plot over second or next deformed images
           [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = PlotstressQuadtree( ...
                DICpara,ResultStrain{ImgSeqNum-1},DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),file_name{1,ImgSeqNum});
 
        end
    end
    
    
    % ------ Save figures for computed stress fields ------
    SaveFigFilesStress; 
    
    % ----- Save strain results ------
    ResultStress{ImgSeqNum-1} = struct('stressxCoord',ResultStrain{ImgSeqNum-1}.strainxCoord,'stressyCoord',ResultStrain{ImgSeqNum-1}.strainyCoord, ...
            'stress_sxx',stress_sxx,'stress_sxy',stress_sxy,'stress_syy',stress_syy, ...
            'stress_principal_max_xyplane',stress_principal_max_xyplane, 'stress_principal_min_xyplane',stress_principal_min_xyplane, ...
            'stress_maxshear_xyplane',stress_maxshear_xyplane,'stress_maxshear_xyz3d',stress_maxshear_xyz3d, ...
            'stress_vonMises',stress_vonMises);
        
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 9 Done ------------ \n \n')

% ------ Save data again including solved stress fields ------
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh',...
                   'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain','ResultStress');





