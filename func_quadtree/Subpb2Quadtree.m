%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function AL-DIC Subproblem 2                             %
% Object: to find deformation field using global methods   %
% Author: Jin Yang                                         %
% Last date modified: 2018.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Uhat] = Subpb2Quadtree(DICmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,winstepsize,waitBarDisplayOrNot)
 
coordinatesFEM = DICmesh.coordinatesFEM;
elementsFEM = DICmesh.elementsFEM;
dirichlet = DICmesh.dirichlet;
neumann = DICmesh.neumann;

DIM = 2; NodesNumPerEle = 4; ClusterNo = 1;
FEMSize = DIM*size(coordinatesFEM,1);
U = USubpb1; F = FSubpb1; W = 0*udual; v = 0*vdual;

%% ====== Initialize variables ======
Uhat = U; U = [U;zeros(DIM*NodesNumPerEle,1)]; v = [v;zeros(DIM*NodesNumPerEle,1)]; 
F = [F;zeros(DIM^2*NodesNumPerEle,1)]; W = [W;zeros(DIM^2*NodesNumPerEle,1)];
% FMinusW1 = F(1:2:end)-W(1:2:end); FMinusW2 = F(2:2:end)-W(2:2:end); % Comment old codes
UMinusv = U-v; FMinusW = F-W;

% ====== Gaussian quadrature parameter ======
% ------ 3*3 Gaussian points ------
gqpt1 = 0; gqpt2 = sqrt(3/5); gqpt3 = -sqrt(3/5); gqpt = [gqpt1,gqpt2,gqpt3];
gqwt1 = 8/9; gqwt2 = 5/9; gqwt3 = 5/9; gqwt = [gqwt1,gqwt2,gqwt3];
% ------ 4*4 Gaussian points ------
% gqpt1 = 0.339981; gqpt2 = -0.339981; gqpt3 = 0.861136; gqpt4 = -0.861136;  
% gqwt1 = 0.652145; gqwt2 = 0.652145; gqwt3 = 0.347855; gqwt4 = 0.347855;  
% gqpt = [gqpt1,gqpt2,gqpt3,gqpt4]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4];
% ------ 5*5 Gaussian points ------
% gqpt1 = 0; gqpt2 = 0.538469; gqpt3 = -0.538469; gqpt4 = 0.90618; gqpt5 = -0.90618;
% gqwt1 = 0.568889; gqwt2 = 0.478629; gqwt3 = 0.478629; gqwt4 = 0.236927; gqwt5 = 0.236927;
% gqpt = [gqpt1,gqpt2,gqpt3,gqpt4,gqpt5]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4,gqwt5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ====== Initialize Global FEM solver ======
ConvergeOrNot = 0; IterStep = 0;
while ConvergeOrNot < 0.5 && IterStep < 1 % Only do one step; IterStep will add value "1" in the first iteration.
    
    IterStep = IterStep+1;
    
    if waitBarDisplayOrNot == 0
        if (ClusterNo==0) || (ClusterNo==1)
            hbar = waitbar(0, ['FE-method to solve Subproblem 2.']);
        else
            hbar=parfor_progressbar(size(elementsFEM,1),['FE-method to solve Subproblem 2.']);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ============= For each element, assemble the stiffness matrix ============
    for eleInd = 1:size(elementsFEM,1) % eleInd is the element index
         
        if waitBarDisplayOrNot == 0
            if ClusterNo== 0 || ClusterNo== 1
                waitbar(eleInd/size(elementsFEM,1));
            else
                hbar.iterate(1);
            end
        end
        
        tempA = zeros(DIM*8,DIM*8); tempb = tempA(:,1);
        
        % ------ Find four corner points ------
        pt1x = coordinatesFEM(elementsFEM(eleInd,1),1); pt1y = coordinatesFEM(elementsFEM(eleInd,1),2);
        pt2x = coordinatesFEM(elementsFEM(eleInd,2),1); pt2y = coordinatesFEM(elementsFEM(eleInd,2),2);
        pt3x = coordinatesFEM(elementsFEM(eleInd,3),1); pt3y = coordinatesFEM(elementsFEM(eleInd,3),2);
        pt4x = coordinatesFEM(elementsFEM(eleInd,4),1); pt4y = coordinatesFEM(elementsFEM(eleInd,4),2);
        
        % ------ Find mid points 5/6/7/8 -------
        if elementsFEM(eleInd,5) ~= 0
            pt5x = coordinatesFEM(elementsFEM(eleInd,5),1); pt5y = coordinatesFEM(elementsFEM(eleInd,5),2);
        else, pt5x = 0; pt5y = 0; 
        end
        if elementsFEM(eleInd,6) ~= 0
            pt6x = coordinatesFEM(elementsFEM(eleInd,6),1); pt6y = coordinatesFEM(elementsFEM(eleInd,6),2);
        else, pt6x = 0; pt6y = 0; 
        end
        if elementsFEM(eleInd,7) ~= 0
            pt7x = coordinatesFEM(elementsFEM(eleInd,7),1); pt7y = coordinatesFEM(elementsFEM(eleInd,7),2);
        else, pt7x = 0; pt7y = 0; 
        end
        if elementsFEM(eleInd,8) ~= 0
            pt8x = coordinatesFEM(elementsFEM(eleInd,8),1); pt8y = coordinatesFEM(elementsFEM(eleInd,8),2);
        else, pt8x = 0; pt8y = 0; 
        end
        
        % ------ Calculate ksi and eta ------
        % lMatrix = [ point1x*point1y point1x point1y 1;
        %             point2x*point2y point2x point2y 1;
        %             point3x*point3y point3x point3y 1;
        %             point4x*point4y point4x point4y 1 ];
          
        % % ------ Find the element nodal indices ------
        % lb = [-1;1;1;-1]; l = linsolve(lMatrix,lb);
        % mb = [-1;-1;1;1]; m = linsolve(lMatrix,mb);
        
        % ------ Find the element nodal indices ------
        tempIndexU = [2*elementsFEM(eleInd,1)-1 2*elementsFEM(eleInd,1) 2*elementsFEM(eleInd,2)-1 2*elementsFEM(eleInd,2) ...
                2*elementsFEM(eleInd,3)-1 2*elementsFEM(eleInd,3) 2*elementsFEM(eleInd,4)-1 2*elementsFEM(eleInd,4) ...
                2*elementsFEM(eleInd,5)-1 2*elementsFEM(eleInd,5) 2*elementsFEM(eleInd,6)-1 2*elementsFEM(eleInd,6) ...
                2*elementsFEM(eleInd,7)-1 2*elementsFEM(eleInd,7) 2*elementsFEM(eleInd,8)-1 2*elementsFEM(eleInd,8)];
        
%         tp = ones(1,DIM);
%         tempIndexU = 2*elementsFEM(j,[tp,2*tp,3*tp,4*tp,5*tp,6*tp,7*tp,8*tp]);
%         tempIndexU(1:2:end) = tempIndexU(1:2:end)-1;  % size of tempIndexU: 1*16

        tempIndexF = [4*elementsFEM(eleInd,1)-3 4*elementsFEM(eleInd,1)-1 4*elementsFEM(eleInd,1)-2 4*elementsFEM(eleInd,1);
                4*elementsFEM(eleInd,1)-3 4*elementsFEM(eleInd,1)-1 4*elementsFEM(eleInd,1)-2 4*elementsFEM(eleInd,1);
                4*elementsFEM(eleInd,2)-3 4*elementsFEM(eleInd,2)-1 4*elementsFEM(eleInd,2)-2 4*elementsFEM(eleInd,2);
                4*elementsFEM(eleInd,2)-3 4*elementsFEM(eleInd,2)-1 4*elementsFEM(eleInd,2)-2 4*elementsFEM(eleInd,2);
                4*elementsFEM(eleInd,3)-3 4*elementsFEM(eleInd,3)-1 4*elementsFEM(eleInd,3)-2 4*elementsFEM(eleInd,3);
                4*elementsFEM(eleInd,3)-3 4*elementsFEM(eleInd,3)-1 4*elementsFEM(eleInd,3)-2 4*elementsFEM(eleInd,3);
                4*elementsFEM(eleInd,4)-3 4*elementsFEM(eleInd,4)-1 4*elementsFEM(eleInd,4)-2 4*elementsFEM(eleInd,4);
                4*elementsFEM(eleInd,4)-3 4*elementsFEM(eleInd,4)-1 4*elementsFEM(eleInd,4)-2 4*elementsFEM(eleInd,4);
                4*elementsFEM(eleInd,5)-3 4*elementsFEM(eleInd,5)-1 4*elementsFEM(eleInd,5)-2 4*elementsFEM(eleInd,5);
                4*elementsFEM(eleInd,5)-3 4*elementsFEM(eleInd,5)-1 4*elementsFEM(eleInd,5)-2 4*elementsFEM(eleInd,5);
                4*elementsFEM(eleInd,6)-3 4*elementsFEM(eleInd,6)-1 4*elementsFEM(eleInd,6)-2 4*elementsFEM(eleInd,6);
                4*elementsFEM(eleInd,6)-3 4*elementsFEM(eleInd,6)-1 4*elementsFEM(eleInd,6)-2 4*elementsFEM(eleInd,6);
                4*elementsFEM(eleInd,7)-3 4*elementsFEM(eleInd,7)-1 4*elementsFEM(eleInd,7)-2 4*elementsFEM(eleInd,7);
                4*elementsFEM(eleInd,7)-3 4*elementsFEM(eleInd,7)-1 4*elementsFEM(eleInd,7)-2 4*elementsFEM(eleInd,7);
                4*elementsFEM(eleInd,8)-3 4*elementsFEM(eleInd,8)-1 4*elementsFEM(eleInd,8)-2 4*elementsFEM(eleInd,8);
                4*elementsFEM(eleInd,8)-3 4*elementsFEM(eleInd,8)-1 4*elementsFEM(eleInd,8)-2 4*elementsFEM(eleInd,8)];
        
        % We don't want temp <= 0, instead, put them to the end
        for tempi = 5:8
            if elementsFEM(eleInd,tempi) == 0
                tempIndexU(2*tempi-1)   = 2*(FEMSize/DIM+1)-1;   tempIndexU(2*tempi)   = 2*(FEMSize/DIM+1);
                tempIndexF(2*tempi-1,1) = 4*(FEMSize/DIM+1)-3;   tempIndexF(2*tempi,1) = 4*(FEMSize/DIM+1)-3;
                tempIndexF(2*tempi-1,2) = 4*(FEMSize/DIM+1)-1;   tempIndexF(2*tempi,2) = 4*(FEMSize/DIM+1)-1;
                tempIndexF(2*tempi-1,3) = 4*(FEMSize/DIM+1)-2;   tempIndexF(2*tempi,3) = 4*(FEMSize/DIM+1)-2;
                tempIndexF(2*tempi-1,4) = 4*(FEMSize/DIM+1);     tempIndexF(2*tempi,4) = 4*(FEMSize/DIM+1);
            end
        end
         
        
        
        % for tempi = 1:size(pointOfx,1)
        %     for tempj = 1:size(pointOfy,1)
        %% ------ Calculate ksi and eta ------
        %        ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
        %        eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
        % %%%%%% Use following codes %%%%%%
        % ------ Nine Gauss integral points ------
        ksietaList = [-sqrt(3/5), 0, sqrt(3/5)]; 
        ksietaWeightList = [5/9, 8/9, 5/9];
        %ksietaList = [-0.90618,-0.538469,0,0.538469,0.90618]; 
        %ksietaWeightList = [0.236927,0.478629,0.568889,0.478629,0.236927];
        
        [ksiAll, etaAll] = ndgrid(ksietaList, ksietaList);
        [weightksiAll, weightetaAll] = ndgrid(ksietaWeightList, ksietaWeightList);
        
        ksiAll = ksiAll(:); etaAll = etaAll(:); weightksiAll = weightksiAll(:); weightetaAll = weightetaAll(:);
        
        
        %for tempi = 1:length(ksietaList)
        %    for tempj = 1:length(ksietaList)
        for tempjj = 1:length(ksiAll)
            
            % ---------------------------
               % ksi = ksietaList(tempi);
               % eta = ksietaList(tempj);
               % weightksi = ksietaWeightList(tempi);
               % weighteta = ksietaWeightList(tempj);
               ksi = ksiAll(tempjj); eta = etaAll(tempjj); weightksi = weightksiAll(tempjj); weighteta = weightetaAll(tempjj);
               
         
                % ------ Calculate N ------
                deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
                if elementsFEM(eleInd,5) ~= 0; deltaPt5 = 1; end
                if elementsFEM(eleInd,6) ~= 0; deltaPt6 = 1; end
                if elementsFEM(eleInd,7) ~= 0; deltaPt7 = 1; end
                if elementsFEM(eleInd,8) ~= 0; deltaPt8 = 1; end
               
                % N5 = deltaPt5*0.5*(1+ksi)*(1-eta^2);
                % N6 = deltaPt6*0.5*(1+eta)*(1-ksi^2);
                % N7 = deltaPt7*0.5*(1-ksi)*(1-eta^2);
                % N8 = deltaPt8*0.5*(1-eta)*(1-ksi^2);
                N5 = deltaPt5*0.5*(1+ksi)*(1-abs(eta));
                N6 = deltaPt6*0.5*(1+eta)*(1-abs(ksi));
                N7 = deltaPt7*0.5*(1-ksi)*(1-abs(eta));
                N8 = deltaPt8*0.5*(1-eta)*(1-abs(ksi));
                
                N1 = (1-ksi)*(1-eta)*0.25 - 0.5*(N7+N8);
                N2 = (1+ksi)*(1-eta)*0.25 - 0.5*(N8+N5);
                N3 = (1+ksi)*(1+eta)*0.25 - 0.5*(N5+N6);
                N4 = (1-ksi)*(1+eta)*0.25 - 0.5*(N6+N7);
                
                
                N = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0;
                         0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8];
                    NDiag = diag([N1,N1,N2,N2,N3,N3,N4,N4,N5,N5,N6,N6,N7,N7,N8,N8]);
                  
                    
                
                % ------ Build J matrix ------
                % Comment: I didn't change Jacobian matrix J when enriched 
                % functions are added.
                J11 = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)*pt1x + funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)*pt2x + ...
                      funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)*pt3x + funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)*pt4x + ...
                      funDN5Dksi(ksi,eta,deltaPt5)*pt5x + funDN6Dksi(ksi,eta,deltaPt6)*pt6x + ...
                      funDN7Dksi(ksi,eta,deltaPt7)*pt7x + funDN8Dksi(ksi,eta,deltaPt8)*pt8x;
                J12 = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)*pt1y + funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)*pt2y + ...
                      funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)*pt3y + funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)*pt4y + ...
                      funDN5Dksi(ksi,eta,deltaPt5)*pt5y + funDN6Dksi(ksi,eta,deltaPt6)*pt6y + ...
                      funDN7Dksi(ksi,eta,deltaPt7)*pt7y + funDN8Dksi(ksi,eta,deltaPt8)*pt8y;
                J21 = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)*pt1x + funDN2Deta(ksi,eta,deltaPt8,deltaPt5)*pt2x + ...
                      funDN3Deta(ksi,eta,deltaPt5,deltaPt6)*pt3x + funDN4Deta(ksi,eta,deltaPt6,deltaPt7)*pt4x + ...
                      funDN5Deta(ksi,eta,deltaPt5)*pt5x + funDN6Deta(ksi,eta,deltaPt6)*pt6x + ...
                      funDN7Deta(ksi,eta,deltaPt7)*pt7x + funDN8Deta(ksi,eta,deltaPt8)*pt8x;
                J22 = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)*pt1y + funDN2Deta(ksi,eta,deltaPt8,deltaPt5)*pt2y + ...
                      funDN3Deta(ksi,eta,deltaPt5,deltaPt6)*pt3y + funDN4Deta(ksi,eta,deltaPt6,deltaPt7)*pt4y + ...
                      funDN5Deta(ksi,eta,deltaPt5)*pt5y + funDN6Deta(ksi,eta,deltaPt6)*pt6y + ...
                      funDN7Deta(ksi,eta,deltaPt7)*pt7y + funDN8Deta(ksi,eta,deltaPt8)*pt8y;
                J = [J11 J12; J21 J22];
                
                Jacobian = det(J);
                InvJ = 1/Jacobian*[J22 -J12; -J21 J11];
                
                
                % ------ Compute DN matrix ------
                DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                        [funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                        funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) 0;
                        funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) 0 ...
                        funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) 0;
                        0 funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                        0 funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                        0 funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) ...
                        0 funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8);
                        0 funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                        0 funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                        0 funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) ...
                        0 funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8)];
                     
                    
                % ------- Comment: Calculate DivFMinusW ---------
                % tempDUDX = DN*FMinusW1(temp);
                % DivFMinusW1 = tempDUDX(1)+tempDUDX(4);
                % tempDUDX = DN*FMinusW2(temp);
                % DivFMinusW2 = tempDUDX(1)+tempDUDX(4);
                % DivFMinusW = [DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2];
                % ------ End of comment of calculating DivFMinusW ------
                
                % ------ Construct A matrix ------
                tempA = tempA + Jacobian*weightksi*weighteta * ( (beta+alpha)*(DN')*(DN) + mu*(N')*N  );
                
                % ------ Construct b vector ------
                tempb = tempb + Jacobian*weightksi*weighteta * ( beta*diag((DN')*(FMinusW(tempIndexF)')) + mu*N'*N*(UMinusv(tempIndexU)) +(alpha)*(DN')*DN*U(tempIndexU) );
                  
            %end
        %end
        end
        
        
        % Store tempA and tempb
        if IterStep==1
            [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
            INDEXAIpar{eleInd}=IndexAXX(:); INDEXAJpar{eleInd}=IndexAYY(:); INDEXAVALpar{eleInd}=tempA(:);
        end
        INDEXBIpar{eleInd}=tempIndexU(:); INDEXBVALpar{eleInd}=tempb(:);   
        
    end
    
    if waitBarDisplayOrNot == 0
      close(hbar);
    end
    
    % ====== Initialize A matrix and b vector ======
    A = sparse(FEMSize+DIM*NodesNumPerEle, FEMSize+DIM*NodesNumPerEle);   b=sparse(FEMSize+DIM*NodesNumPerEle,1);
    for eleInd = 1:size(elementsFEM,1)
        A = A + sparse(INDEXAIpar{eleInd}, INDEXAJpar{eleInd}, INDEXAVALpar{eleInd}, FEMSize+DIM*NodesNumPerEle , FEMSize+DIM*NodesNumPerEle) ;
        b = b + sparse(INDEXBIpar{eleInd},ones(length(INDEXBIpar{eleInd}),1),INDEXBVALpar{eleInd}, FEMSize+DIM*NodesNumPerEle , 1) ;
    end
% "DIM*NodesNumPerEle" is because I put all the zeros in elementsFEM to the end, and in NC function there are "DIM*NodesNumPerEle" entries




    %if IterStep == 1
        % AMatrixRegularized = A + 1e-3*max(diag(A))*speye(size(A,1),size(A,2));
        % AMatrixRegularized = A; tempVal = 1e-3*max(diag(A));
        % for tempi = 1:size(A,1)
        %     AMatrixRegularized(tempi,tempi) = AMatrixRegularized(tempi,tempi) + tempVal;
        % send
    %end
    
    coordsIndexInvolved = unique(elementsFEM(:,1:8)); % Need modification for triangle elementsFEM
     
    UIndexInvolved = zeros(2*(length(coordsIndexInvolved)-1),1);
    % Not including the first 0-th entry
    for tempi = 1:(size(coordsIndexInvolved,1)-1)
        UIndexInvolved(2*tempi-1:2*tempi) = [2*coordsIndexInvolved(tempi+1)-1; 2*coordsIndexInvolved(tempi+1)];
    end
    
    % ========= Set Dirichlet and Neumann boundary conditions =========
    if isempty(dirichlet) ~= 1
        dirichlettemp = [2*dirichlet(:); 2*dirichlet(:)-1];
    else
        dirichlettemp = [];
    end
%     if isempty(neumann) ~= 1
%         neumanntemp = [2*neumann(:,1); 2*neumann(:,1)-1; 2*neumann(:,2); 2*neumann(:,2)-1];
%     else
%         neumanntemp = [];
%     end
    FreeNodes = setdiff(UIndexInvolved,unique([dirichlettemp]));
    
    % ========= Neumann conditions ===========
    % Last step boundary condition force
    BCForce = - 1/winstepsize * FSubpb1;
    for tempj = 1:size(neumann,1)
        b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + 0.5*norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
          *( ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) ) );
        b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + 0.5*norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
         * (( BCForce(4*neumann(tempj,1:2)-2) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) ) );
    end
    
    % ========= Dirichlet conditions ==========
    UhatOld = Uhat; Uhat = sparse(FEMSize+DIM*NodesNumPerEle,1);
    for tempi = 1:DIM
        Uhat(DIM*unique(dirichlet)-(tempi-1)) = U(DIM*unique(dirichlet)-(tempi-1)); 
    end
    b = b - A * Uhat;
     
    % ========= Solve FEM problem ===========
    Uhat(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);
    UhatNew = Uhat(1:end-DIM*NodesNumPerEle);
     
    if norm(UhatNew-UhatOld)/sqrt(length(UhatOld)) < 1e-3
        ConvergeOrNot = 1;
    end
    
    Uhat = UhatNew;
    
    
    
end

end


%% ========= subroutines for  FEM Q4 shape function derivatives ========
function DN5Dksi=funDN5Dksi(ksi,eta,deltaPt5)
DN5Dksi = deltaPt5*0.5*(1-abs(eta)) ;
%DN5Dksi = deltaPt5*0.5*(1-eta^2);
end
function DN5Deta=funDN5Deta(ksi,eta,deltaPt5)
DN5Deta = deltaPt5*0.5*(1+ksi)*sign(-eta);
%DN5Deta = deltaPt5*(-1)*(1+ksi)*eta;
end
function DN6Dksi=funDN6Dksi(ksi,eta,deltaPt6)
DN6Dksi = deltaPt6*0.5*(1+eta)*sign(-ksi);
%DN6Dksi = deltaPt6*(-1)*(1+eta)*ksi;
end
function DN6Deta=funDN6Deta(ksi,eta,deltaPt6)
DN6Deta = deltaPt6*0.5*(1-abs(ksi));
%DN6Deta = deltaPt6*0.5*(1-ksi^2);
end
function DN7Dksi=funDN7Dksi(ksi,eta,deltaPt7)
DN7Dksi = deltaPt7*0.5*(-1)*(1-abs(eta));
%DN7Dksi = deltaPt7*(-0.5)*(1-eta^2);
end
function DN7Deta=funDN7Deta(ksi,eta,deltaPt7)
DN7Deta = deltaPt7*0.5*(1-ksi)*sign(-eta);
%DN7Deta = deltaPt7*(-1)*(1-ksi)*eta;
end
function DN8Dksi=funDN8Dksi(ksi,eta,deltaPt8)
DN8Dksi = deltaPt8*0.5*(1-eta)*sign(-ksi);
%DN8Dksi = deltaPt8*(-1)*(1-eta)*ksi;
end
function DN8Deta=funDN8Deta(ksi,eta,deltaPt8)
DN8Deta = deltaPt8*0.5*(-1)*(1-abs(ksi));
%DN8Deta = deltaPt8*(-0.5)*(1-ksi^2);
end
function DN1Dksi = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)
DN1Dksi = -0.25*(1-eta)-0.5*( funDN7Dksi(ksi,eta,deltaPt7) + funDN8Dksi(ksi,eta,deltaPt8) );
end
function DN1Deta = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)
DN1Deta = -0.25*(1-ksi)-0.5*( funDN7Deta(ksi,eta,deltaPt7) + funDN8Deta(ksi,eta,deltaPt8) );
end
function DN2Dksi = funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)
DN2Dksi = 0.25*(1-eta)-0.5*( funDN8Dksi(ksi,eta,deltaPt8) + funDN5Dksi(ksi,eta,deltaPt5) );
end
function DN2Deta = funDN2Deta(ksi,eta,deltaPt8,deltaPt5)
DN2Deta = -0.25*(1+ksi)-0.5*( funDN8Deta(ksi,eta,deltaPt8) + funDN5Deta(ksi,eta,deltaPt5) );
end
function DN3Dksi = funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)
DN3Dksi = 0.25*(1+eta)-0.5*( funDN5Dksi(ksi,eta,deltaPt5) + funDN6Dksi(ksi,eta,deltaPt6) );
end
function DN3Deta = funDN3Deta(ksi,eta,deltaPt5,deltaPt6)
DN3Deta = 0.25*(1+ksi)-0.5*( funDN5Deta(ksi,eta,deltaPt5) + funDN6Deta(ksi,eta,deltaPt6) );
end
function DN4Dksi = funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)
DN4Dksi = -0.25*(1+eta)-0.5*( funDN6Dksi(ksi,eta,deltaPt6) + funDN7Dksi(ksi,eta,deltaPt7) );
end
function DN4Deta = funDN4Deta(ksi,eta,deltaPt6,deltaPt7)
DN4Deta = 0.25*(1-ksi)-0.5*( funDN6Deta(ksi,eta,deltaPt6) + funDN7Deta(ksi,eta,deltaPt7) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 


