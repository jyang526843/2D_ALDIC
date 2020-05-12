% ---------------------------------------------
% Augmented Lagrangian DIC
% Author: Jin Yang, PhD. Caltech
% Contact: yangjin@caltech.edu
% 2015.04,06,07; 2016.03,04
% ---------------------------------------------

%% Section 1: Clear MATLAB environment & mex set up Spline interpolation  
close all; clear; clc; 
fprintf('------------ Section 1 Start ------------ \n')
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
% % cd("./Splines_interp/lib_matlab"); CompileLib; cd("../../");  % % mex bi-cubic spline interpolations
% % addpath("./Splines_interp/lib_matlab"); % dbstop if error % % Old version codes.
mex -O ba_interp2.cpp; 
addpath("./func"); addpath("./src"); addpath("./plotFiles/"); addpath("./plotFiles/export_fig-d966721/");
% addpath("./YOUR IMAGE FOLDER"); 
fprintf('------------ Section 1 Done ------------ \n \n')


%% Section 2: Load DIC parameters and set up DIC parameters 
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ======
[file_name,Img,DICpara] = ReadImage; close all;
% %%%%%% Uncomment the line below to change the DIC computing ROI manually %%%%%%
%gridxROIRange = [gridxROIRange1,gridxROIRange2]; gridyROIRange = [Val1, Val2];
%gridxROIRange = [224,918]; gridyROIRange = [787,1162];
% ====== Normalize images ======
[ImgNormalized,DICpara.gridxyROIRange] = funNormalizeImg(Img,DICpara.gridxyROIRange);  
% ====== Initialize variable storage ======
ResultDisp = cell(length(ImgNormalized)-1,1); 
ResultDefGrad = cell(length(ImgNormalized)-1,1);
ResultStrain = cell(length(ImgNormalized)-1,1);
ResultFEMesh = cell(ceil((length(ImgNormalized)-1)/DICpara.ImgSeqIncUnit),1); % For incremental DIC mode
fprintf('------------ Section 2 Done ------------ \n \n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fNormalized = ImgNormalized{1}; 
for ImgSeqNum = 2:length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %% Section 3:
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update initial guess for ALDIC
    % The key idea is to either to use FFT peak fitting, or use last frame
    % results for the next new frame;
    % Particularly in incremental mode, the reference image can also be updated.
    % fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)};
    gNormalized = ImgNormalized{ImgSeqNum}; NewFFTSearchCheck = 0; DICpara.NewFFTSearch = 0;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while NewFFTSearchCheck == 0 
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1
        % ====== Integer Search ======
        % Old version: 
        [DICpara.SizeOfFFTSearchRegion,x0temp,y0temp,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara);
        % New version: adaptive search initial guess
        % [x0temp,y0temp,u,v,cc]= IntegerSearchMg(fNormalized,gNormalized,file_name,DICpara);
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

    
    %% Section 4
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
    close all; USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); 
    Plotuv(USubpb1World,DICmesh.x0,DICmesh.y0World);
    disp('--- Start to manually remove bad points ---')
    u = reshape(USubpb1(1:2:end),size(DICmesh.x0,1),size(DICmesh.x0,2)); v = reshape(USubpb1(2:2:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    [u,v,~,Local_BadptRow,Local_BadptCol,RemoveOutliersList] = funRemoveOutliers(u',v',[],0.5,100); u=u';v=v';
    USubpb1(1:2:end) = reshape(u,size(DICmesh.coordinatesFEM,1),1); USubpb1(2:2:end) = reshape(v,size(DICmesh.coordinatesFEM,1),1);
    f11 = reshape(FSubpb1(1:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2)); f21 = reshape(FSubpb1(2:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    f12 = reshape(FSubpb1(3:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2)); f22 = reshape(FSubpb1(4:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
    f11=f11'; f21=f21'; f12=f12'; f22=f22';
    f11(RemoveOutliersList) = NaN; f21(RemoveOutliersList) = NaN; f12(RemoveOutliersList) = NaN; f22(RemoveOutliersList) = NaN;
    f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4); f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
    f11=f11'; f21=f21'; f12=f12'; f22=f22';
    FSubpb1(1:4:end) = reshape(f11,size(DICmesh.coordinatesFEM,1),1); FSubpb1(2:4:end) = reshape(f21,size(DICmesh.coordinatesFEM,1),1);
    FSubpb1(3:4:end) = reshape(f12,size(DICmesh.coordinatesFEM,1),1); FSubpb1(4:4:end) = reshape(f22,size(DICmesh.coordinatesFEM,1),1);
    disp('--- Remove bad points done ---')
    % Comment END
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % You won't use these codes now. Please ignore all the codes related
    % with NewFFTSearchCheck and don't change them.
    % --- --- Check: need redo FFT search or not ------
    %if (LocalICGNBadPtNumtemp/size(DICmesh.coordinatesFEM,1) > 0.1) %&& (mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0) % 10% bad points
    %    DICpara.NewFFTSearch = 1;
    %else
        NewFFTSearchCheck = 1; %DICpara.NewFFTSearch = 0;
    %end
    end % Please ignore this while-loop now.   
    
    % ------ Plot ------
    USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); 
    FSubpb1World = FSubpb1; % FSubpb1World(2:2:end) = -FSubpb1(2:2:end);
    %close all; Plotuv(USubpb1World,x0,y0World); Plotdisp_show(USubpb1World,coordinatesFEMWorld,elementsFEM);
    %Plotstrain_show(FSubpb1World,coordinatesFEMWorld,elementsFEM);
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    fprintf('------------ Section 4 Done ------------ \n \n')

   
    %% Section 5
    fprintf('------------ Section 5 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve global step in ALDIC: Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; LevelNo=1;
    FSubpb1 = funSmoothStrain(FSubpb1,DICmesh,DICpara);
    
	% ====== Define penalty parameter ======
    mu = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1; 
    betaList = [1e-3,sqrt(1e-5),1e-2,sqrt(1e-3),1e-1,sqrt(1e-1)]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
    % ====== Check to use FD or FE methods to solve Subpb2 step ======
	if DICpara.Subpb2FDOrFEM == 1 % Using FD method
		% ====== Build sparse finite difference operator ======
		disp('Assemble finite difference operator D');
		M = size(DICmesh.x0,1); N = size(DICmesh.x0,2);
		tic; Rad = 1; D = funDerivativeOp((M-2*Rad),(N-2*Rad),DICpara.winstepsize); % D = sparse(4*(M-2*Rad)*(N-2*Rad), 2*(M-2*Rad)*(N-2*Rad));
		D2 = funDerivativeOp(M,N,DICpara.winstepsize); toc
		disp('Finish assembling finite difference operator D');
		% ===== Solver using finite difference approximation ======
		tic; a = FSubpb1-udual; b = USubpb1-vdual; 
		Rad = 1; [temp3,temp4] = funFDNeumannBCInd(size(DICmesh.coordinatesFEM,1),M,N,Rad); % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
		atemp = a(temp3); btemp = b(temp4); hbar = waitbar(0,'Please wait for Subproblem 2 global step!');  
        for tempk = 1:length(betaList)
            beta = betaList(tempk);
            tempAMatrixSub2 = (beta*(D')*D) + mu*speye(2*(M-2*Rad)*(N-2*Rad));
            USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp ) ;
            USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp; 
            FSubpb2 = D2*USubpb2; 
            
            Err1(tempk) = norm(USubpb1-USubpb2,2);
            Err2(tempk) = norm(FSubpb1-FSubpb2,2);
            waitbar(tempk/(length(betaList)+1));
        end
        ErrSum = Err1+Err2*mean(DICpara.winstepsize)^2; [~,indexOfbeta] = min(ErrSum);
        try
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch
            beta = betaList(indexOfbeta);
        end
        % Using optimal beta to solve again
        tempAMatrixSub2 = (beta*(D')*D) + mu*speye(2*(M-2*Rad)*(N-2*Rad));
        USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp) ;
        USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp;
        waitbar(1); close(hbar);
		%%%%%%%%%%%%%% End of using finite difference approximation %%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
	else %Subpb2FDOrFEM: Using FE method
        M = size(DICmesh.x0,1); N = size(DICmesh.x0,2); GaussPtOrder = 2; alpha = 0;
        close all; hbar = waitbar(0,'Please wait for Subproblem 2 global step!');  
        % ====== Solver using finite element method ======
        for tempk = 1:length(betaList)
            beta = betaList(tempk);
            GaussPtOrder = 2; alpha = 0; [USubpb2] = Subpb2(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder);
            [FSubpb2] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,USubpb2,GaussPtOrder);
            Err1(tempk) = norm(USubpb1-USubpb2,2);
            Err2(tempk) = norm(FSubpb1-FSubpb2,2);
            waitbar(tempk/(length(betaList)+1));
        end
        Err1Norm = (Err1-mean(Err1))/std(Err1); figure, plot(Err1Norm);
        Err2Norm = (Err2-mean(Err2))/std(Err2); figure, plot(Err2Norm);
        ErrSum = Err1Norm+Err2Norm; figure,plot(ErrSum); [~,indexOfbeta] = min(ErrSum);
        try
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch
            beta = betaList(indexOfbeta);
        end
        % Using optimal beta to solve again
        [USubpb2] = Subpb2(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder);
        USubpb2 = full(USubpb2);
        waitbar(1); close(hbar);
    end
	ALSub2Time(ALSolveStep) = toc; toc
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------- Before computing strain, we smooth the displacement field -------
    %USubpb2 = funSmoothDisp(USubpb2,coordinatesFEM,elementsFEM,x0,y0,winstepsize,DispFilterSize,DispFilterStd);
    % ------- Compute strain field --------
    if DICpara.Subpb2FDOrFEM == 1 %FD 
		FSubpb2 = D2*USubpb2; % D2 = funDerivativeOp(M,N,winstepsize);
	else %FEM
		[FSubpb2] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,USubpb2,GaussPtOrder);
    end
    
    % ------- Smooth strain field --------
    FSubpb2 = funSmoothStrain(FSubpb2,DICmesh,DICpara);

    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end); FSubpb2World = FSubpb2;  
    close all; Plotuv(USubpb2World,DICmesh.x0,DICmesh.y0World); Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);
    %Plotstrain_show(FSubpb2World,coordinatesFEMWorld,elementsFEM);

    % ======= Update dual variables =======
    if DICpara.Subpb2FDOrFEM == 1 %FD
		udualtemp1 = (FSubpb2 - FSubpb1); udualtemp2 = udualtemp1(temp3);
		vdualtemp1 = (USubpb2 - USubpb1); vdualtemp2 = vdualtemp1(temp4);
		udual = zeros(4*M*N,1); vdual = zeros(2*M*N,1);
		udual(temp3) = udualtemp2; vdual(temp4) = vdualtemp2;
	else  % FEM or other methods
		udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
	end
    save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Section 6
    fprintf('------------ Section 6 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to run ADMM iteration: Subproblem 1 & 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ==================== Global AL Loop ==========================
    ALSolveStep = 1; tol2 = 1e-4; UpdateY = 1e4; CrackOrNot = 0; CrackPath1 = [0,0]; CrackPath2 = [0,0]; CrackTip = [0,0]; 
    HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end

    while (ALSolveStep < 5)
        ALSolveStep = ALSolveStep + 1;  % Update using the last step
        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
        tic;[USubpb1,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1(USubpb2,FSubpb2,udual,vdual,DICmesh.coordinatesFEM,...
        Df,fNormalized,gNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        FSubpb1 = FSubpb2; toc 
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp;  ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------  Manually find some bad points from Local Subset ICGN step ------
        % disp('--- Start to manually remove bad points --- \n')
        % disp('    Comment codes here if you do not have bad local points \n')
        % % Comment START
        % close all; Plotuv(USubpb1,x0,y0);
        % u = reshape(USubpb1(1:2:end),M,N); v = reshape(USubpb1(2:2:end),M,N);
        % [u,v,~,Subpb1_BadptRow,Subpb1_BadptCol] = funRemoveOutliers(u',v',[],0.5,100,Local_BadptRow,Local_BadptCol); u=u';v=v';
        % disp('--- Remove bad points done ---')
        % USubpb1(1:2:end) = reshape(u,size(elements,1),1); USubpb1(2:2:end) = reshape(v,size(elements,1),1);
        % close all; Plotuv(USubpb1,x0,y0); Plotdisp_show(USubpb1,coordinatesFEM,elementsFEM);
        % Comment END
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
        USubpb1 = funSmoothDisp(USubpb1,DICmesh,DICpara);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ============== Subproblem 2 ==============
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        if DICpara.Subpb2FDOrFEM == 1 %FD
			% ------- using finite difference approximation --------
			tic; a = FSubpb1-udual; b = USubpb1-vdual; atemp = a(temp3); btemp = b(temp4);
			USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp ) ;
            USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp; %toc
			% ------- End of using finite difference approximation --------
		else % FEM 
			tic; [USubpb2] = Subpb2(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder);
            USubpb2 = full(USubpb2);
		end
		ALSub2Time(ALSolveStep) = toc; toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ------- Before computing strain, we smooth the displacement field -------
        USubpb2 = funSmoothDisp(USubpb2,DICmesh,DICpara);
        % ------- Compute strain field --------
        if DICpara.Subpb2FDOrFEM == 1 %FD
			FSubpb2 = D2*USubpb2; % D2 = funDerivativeOp(M,N,winstepsize);
		else %FEM
			GaussPtOrder = 2; [FSubpb2] = funGlobal_NodalStrainAvg(DICmesh.coordinatesFEM,DICmesh.elementsFEM,USubpb2,GaussPtOrder);
		end
		
		% ------- Smooth strain field --------
        FSubpb2 = funSmoothStrain(FSubpb2,DICmesh,DICpara);
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
        if DICpara.Subpb2FDOrFEM == 1 %FD
			udualtemp1 =  (FSubpb2 - FSubpb1); udualtemp2 = udualtemp1(temp3);
			vdualtemp1 =  (USubpb2 - USubpb1); vdualtemp2 = vdualtemp1(temp4);
			udual(temp3) = udual(temp3)+udualtemp2; 
			vdual(temp4) = vdual(temp4)+vdualtemp2;
		else %FEM
			udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1; 
		end

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


% ------ Plot ------
USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
FSubpb2World = FSubpb2; % FSubpb1World(2:2:end) = -FSubpb1(2:2:end);
close all; Plotuv(USubpb2World,DICmesh.x0,DICmesh.y0World); 
Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM);
% Plotstrain_show(FSubpb2World,coordinatesFEMWorld,elementsFEM);

%% ------ Save results ------
% Find img name and save all the results 
[~,imgname,imgext] = fileparts(file_name{1,end});
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh');

  
%% Section 7
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

%% Section 8
fprintf('------------ Section 8 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain and plot figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Smooth displacements ------
DICpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');
% ------ Choose strain computation method ------
DICpara.MethodToComputeStrain = funParaInput('StrainMethodOp'); 
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DICpara.StrainType = funParaInput('StrainType');
% ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat');
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xList = min(coordinatesFEM(:,1)):DICpara.winstepsize:max(coordinatesFEM(:,1)); M = length(xList);
    yList = min(coordinatesFEM(:,2)):DICpara.winstepsize:max(coordinatesFEM(:,2)); N = length(yList);
    %[x0,y0] = ndgrid(xList,yList); 
    [x0,y0] = ndgrid(xList,yList); 
    x0 = x0-reshape(UFEMesh(1:2:end),size(x0,1),size(x0,2)); 
    y0 = y0-reshape(UFEMesh(2:2:end),size(y0,1),size(y0,2)); 
    y0World = (size(ImgNormalized{1},2)+1-y0); 
    coordinatesFEMWorld = [coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)];

    % ------ Plotting and Compute Strain-------
    M = size(x0,1); N = size(x0,2);
    if size(USubpb2,1) == 1
        ULocal = USubpb2_New.USubpb2; FLocal = FSubpb2.FSubpb2; 
    else
        ULocal = USubpb2; FLocal = FSubpb2;
    end
    % ULocal(1:2:end)= -ULocal(1:2:end); 
    % FLocal(1:4:end)= -1*FLocal(1:4:end); FLocal(3:4:end)= -1*FLocal(3:4:end); % because u is flipped sign
    UWorld = ULocal; UWorld(2:2:end) = -UWorld(2:2:end); FWorld = FLocal; close all; Plotuv(UWorld,x0,y0World);
    % tic; D = funDerivativeOp(M,N,winstepsize); toc

    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    DispFilterSize=0; DispFilterStd=0; SmoothTimes = 0;
    try
        while DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            ULocal = funSmoothDisp(ULocal,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
    
    % ----- Compute strain field ------
    ComputeStrain; % Compute strain
    % % ------- Add filter and plot strain field -------
    %Plotstrain_Fij;

    % ------ Plot disp and strain ------
    if DICpara.OrigDICImgTransparency == 0
        Plotdisp_show(UWorld,coordinatesFEMWorld,elementsFEM);
        Plotstrain0(FStraintemp,x0(1+Rad:M-Rad,1+Rad:N-Rad),y0(1+Rad:M-Rad,1+Rad:N-Rad),size(ImgNormalized{1}),...
        file_name{1,ImgSeqNum},DICpara.OrigDICImgTransparency); 
    else
        Plotdisp(UWorld,x0,y0,size(ImgNormalized{1}),file_name{1,ImgSeqNum},DICpara.OrigDICImgTransparency);
        Plotstrain(UWorld,Rad,FStraintemp,x0(1+Rad:M-Rad,1+Rad:N-Rad),y0(1+Rad:M-Rad,1+Rad:N-Rad),size(ImgNormalized{1}),...
        file_name{1,ImgSeqNum},DICpara.OrigDICImgTransparency); 
    end

    ResultStrain{ImgSeqNum-1}.Strain = FStraintemp;
    
    % ------ Save figures ------
    SaveFigFiles;
    
end
% ------ Save movies ------
fprintf('------------ Section 8 Done ------------ \n \n')

%% ------- Store images ------
% imgNo = 0570;
% mean(FStraintemp(1:4:end)) % -0.0029(57939_WS20ST10) -0.0029(57939_WS40ST10) -5.03e-4(57383_WS20ST10)  -5.02e-4(57383_WS40ST10)
% mean(FStraintemp(4:4:end)) % 0.0056(57939_WS20ST10) 0.0056(57939_WS40ST10) 3.804e-4(57383_WS20ST10)  3.85e-4(57383_WS20ST10)
% 0.5*(mean(FStraintemp(2:4:end))+mean(FStraintemp(3:4:end))) % -5.7438e-4(57939_WS20ST10)  -5.7665e-4(57939_WS40ST10) -1.06e-4(57383_WS20ST10)  -1.0456e-4(57383_WS40ST10)

%% Save data again including strain solve method
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultStrain','ResultFEMesh',...
    'ALSub1Time','ALSub2Time','ALSolveStep');





