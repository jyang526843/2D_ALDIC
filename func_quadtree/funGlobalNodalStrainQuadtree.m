function [Strain,StrainGaussPt,CoordGaussPt] = funGlobalNodalStrainQuadtree(DICmesh,U,GaussPtOrder,waitBarDisplayOrNot)

coordinatesFEM = DICmesh.coordinatesFEM; 
elementsFEM = DICmesh.elementsFEM;
h = DICmesh.elementMinSize;

DIM = 2; NodesNumPerEle = 4;
SizeCoords = size(coordinatesFEM,1);

%% ====== Initialize strain Gauss points ======
Strain = 0; U = [U;zeros(DIM*NodesNumPerEle,1)]; % For zero entries in elements
FStrainAvgTimes = zeros(4*SizeCoords,1); FStrain = zeros(4*SizeCoords,1);
StrainGaussPt = zeros(4*size(elementsFEM,1),4); CoordGaussPt = zeros(4*size(elementsFEM,1),2);

% ====== Gaussian quadrature parameter ======
switch GaussPtOrder
    case 2 % ------ 2*2 Gauss points ------
        gqpt1 = -1/sqrt(3); gqpt2 = 1/sqrt(3); gqpt = [gqpt1,gqpt2]; 
        gqwt1 = 1; gqwt2 = 1; gqwt = [gqwt1,gqwt2];
    case 3 % ------ 3*3 Gauss points ------
        gqpt1 = 0; gqpt2 = sqrt(3/5); gqpt3 = -sqrt(3/5); gqpt = [gqpt1,gqpt2,gqpt3];
        gqwt1 = 8/9; gqwt2 = 5/9; gqwt3 = 5/9; gqwt = [gqwt1,gqwt2,gqwt3];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waitBarDisplayOrNot == 0
   hbar = waitbar(0, ['Update deformation gradient tensor F after solving Subproblem 2.']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for eleInd = 1:size(elementsFEM,1) % eleInd is the element index
    
    if waitBarDisplayOrNot == 0
    waitbar(eleInd/size(elementsFEM,1));
    end
    
    StrainWithinEachElementGausspoint = zeros(4,4);
    
    % ----- Find four corner points -------
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

    % ------ Calculate ksi and eta --------
    lMatrix = [ pt1x*pt1y pt1x pt1y 1;
                pt2x*pt2y pt2x pt2y 1;
                pt3x*pt3y pt3x pt3y 1;
                pt4x*pt4y pt4x pt4y 1 ];

    % ------ Find the element nodal indices ------
    lb = [-1;1;1;-1]; l = linsolve(lMatrix,lb);
    mb = [-1;-1;1;1]; m = linsolve(lMatrix,mb);

    % ------ Find the element nodal indices ------
    tempIndexU = [2*elementsFEM(eleInd,1)-1 2*elementsFEM(eleInd,1) 2*elementsFEM(eleInd,2)-1 2*elementsFEM(eleInd,2) ...
            2*elementsFEM(eleInd,3)-1 2*elementsFEM(eleInd,3) 2*elementsFEM(eleInd,4)-1 2*elementsFEM(eleInd,4) ...
            2*elementsFEM(eleInd,5)-1 2*elementsFEM(eleInd,5) 2*elementsFEM(eleInd,6)-1 2*elementsFEM(eleInd,6) ...
            2*elementsFEM(eleInd,7)-1 2*elementsFEM(eleInd,7) 2*elementsFEM(eleInd,8)-1 2*elementsFEM(eleInd,8)];
        
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elementsFEM(eleInd,tempi) == 0
            tempIndexU(2*tempi-1) = 2*(SizeCoords+1)-1;    
            tempIndexU(2*tempi)   = 2*(SizeCoords+1);
        end
    end
     
    % ------ Set Gauss points ------
    pt1 = elementsFEM(eleInd,1); pt2 = elementsFEM(eleInd,3);
    ptOfx = zeros(length(gqwt),1); ptOfy = zeros(length(gqwt),1);
    for tempi = 1:length(gqwt)
       ptOfx(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,1)-coordinatesFEM(pt1,1))+0.5*(coordinatesFEM(pt2,1)+coordinatesFEM(pt1,1));
       ptOfy(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,2)-coordinatesFEM(pt1,2))+0.5*(coordinatesFEM(pt2,2)+coordinatesFEM(pt1,2));
    end
    for tempi = 1:length(ptOfx)
       for tempj = 1:length(ptOfy)
    % ------ Start four Gauss points ------
    %for tempi = 1:2
    %    for tempj = 1:2
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*ptOfx(tempi)*ptOfy(tempj) + l(2)*ptOfx(tempi) + l(3)*ptOfy(tempj) + l(4) ;
            eta = m(1)*ptOfx(tempi)*ptOfy(tempj) + m(2)*ptOfx(tempi) + m(3)*ptOfy(tempj) + m(4) ;
            % ksi = 2*tempj-3; eta = 2*tempi-3;
            % if (tempi == 1) && (tempj == 1)
            %     ksi = -1/sqrt(3); eta = -1/sqrt(3);
            % elseif (tempi == 1) && (tempj == 2)
            %     ksi = 1/sqrt(3); eta = -1/sqrt(3);
            % elseif (tempi == 2) && (tempj == 1)
            %     ksi = 1/sqrt(3); eta = 1/sqrt(3);
            % elseif (tempi == 2) && (tempj == 2)
            %     ksi = -1/sqrt(3); eta = 1/sqrt(3);
            % end
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elementsFEM(eleInd,5) ~= 0; deltaPt5 = 1; end
            if elementsFEM(eleInd,6) ~= 0; deltaPt6 = 1; end
            if elementsFEM(eleInd,7) ~= 0; deltaPt7 = 1; end
            if elementsFEM(eleInd,8) ~= 0; deltaPt8 = 1; end

            N5 = deltaPt5*0.5*(1+ksi)*(1-abs(eta));
            N6 = deltaPt6*0.5*(1+eta)*(1-abs(ksi));
            N7 = deltaPt7*0.5*(1-ksi)*(1-abs(eta));
            N8 = deltaPt8*0.5*(1-eta)*(1-abs(ksi));

            N1 = (1-ksi)*(1-eta)*0.25 - 0.5*(N7+N8);
            N2 = (1+ksi)*(1-eta)*0.25 - 0.5*(N8+N5);
            N3 = (1+ksi)*(1+eta)*0.25 - 0.5*(N5+N6);
            N4 = (1-ksi)*(1+eta)*0.25 - 0.5*(N6+N7);

             
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
                    0 funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8);];
             
            StrainWithinEachElementGausspoint(length(gqpt)*(tempi-1)+tempj,1:4) = DN*U(tempIndexU);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Replace with 9-Gauss point %%%%%%%%%%%%%%%%%%%%%%
    switch GaussPtOrder
        case 3
            StrainGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:4) = StrainWithinEachElementGausspoint;
            CoordGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:2) = ...
                    [ptOfx(1),ptOfy(1);ptOfx(1),ptOfy(2);ptOfx(1),ptOfy(3);
                     ptOfx(2),ptOfy(1);ptOfx(2),ptOfy(2);ptOfx(2),ptOfy(3);
                     ptOfx(3),ptOfy(1);ptOfx(3),ptOfy(2);ptOfx(3),ptOfy(3)];
        %%%%%%%%%%%%%%%%%%% Comment following 4-Gauss point %%%%%%%%%%%%%%%%%%%%%%
        case 2
            StrainGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:4) = StrainWithinEachElementGausspoint;
            CoordGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:2) = ...
                    [ptOfx(1),ptOfy(1);ptOfx(1),ptOfy(2); 
                     ptOfx(2),ptOfy(1);ptOfx(2),ptOfy(2)];
                 
%             MatrixExtrapolation = [1+0.5*sqrt(3)  -0.5           1-0.5*sqrt(3)   -0.5;
%                                     -0.5           1+0.5*sqrt(3)  -0.5            1-0.5*sqrt(3);
%                                     1-0.5*sqrt(3)  -0.5           1+0.5*sqrt(3)   -0.5;
%                                     -0.5           1-0.5*sqrt(3)  -0.5            1+0.5*sqrt(3)];
%     
%             % ------ Nodal points strain extrapolation using Gauss points -----
%             StrainExxWithinEachElementNodalpoint =  MatrixExtrapolation * ...
%                 [StrainWithinEachElementGausspoint(1,1);
%                 StrainWithinEachElementGausspoint(2,1);
%                 StrainWithinEachElementGausspoint(3,1);
%                 StrainWithinEachElementGausspoint(4,1)];
% 
%             StrainExyWithinEachElementNodalpoint = MatrixExtrapolation * ...
%                 [StrainWithinEachElementGausspoint(1,2);
%                 StrainWithinEachElementGausspoint(2,2);
%                 StrainWithinEachElementGausspoint(3,2);
%                 StrainWithinEachElementGausspoint(4,2)];
% 
%             StrainEyxWithinEachElementNodalpoint = MatrixExtrapolation * ...
%                 [StrainWithinEachElementGausspoint(1,3);
%                 StrainWithinEachElementGausspoint(2,3);
%                 StrainWithinEachElementGausspoint(3,3);
%                 StrainWithinEachElementGausspoint(4,3)];
% 
%             StrainEyyWithinEachElementNodalpoint = MatrixExtrapolation * ...
%                 [StrainWithinEachElementGausspoint(1,4);
%                 StrainWithinEachElementGausspoint(2,4);
%                 StrainWithinEachElementGausspoint(3,4);
%                 StrainWithinEachElementGausspoint(4,4)];
% 
%             StrainWithinEachElementGausspoint(1,1) = StrainExxWithinEachElementNodalpoint(1);
%             StrainWithinEachElementGausspoint(2,1) = StrainExxWithinEachElementNodalpoint(2);
%             StrainWithinEachElementGausspoint(3,1) = StrainExxWithinEachElementNodalpoint(3);
%             StrainWithinEachElementGausspoint(4,1) = StrainExxWithinEachElementNodalpoint(4);
% 
%             StrainWithinEachElementGausspoint(1,2) = StrainExyWithinEachElementNodalpoint(1);
%             StrainWithinEachElementGausspoint(2,2) = StrainExyWithinEachElementNodalpoint(2);
%             StrainWithinEachElementGausspoint(3,2) = StrainExyWithinEachElementNodalpoint(3);
%             StrainWithinEachElementGausspoint(4,2) = StrainExyWithinEachElementNodalpoint(4);
% 
%             StrainWithinEachElementGausspoint(1,3) = StrainEyxWithinEachElementNodalpoint(1);
%             StrainWithinEachElementGausspoint(2,3) = StrainEyxWithinEachElementNodalpoint(2);
%             StrainWithinEachElementGausspoint(3,3) = StrainEyxWithinEachElementNodalpoint(3);
%             StrainWithinEachElementGausspoint(4,3) = StrainEyxWithinEachElementNodalpoint(4);
% 
%             StrainWithinEachElementGausspoint(1,4) = StrainEyyWithinEachElementNodalpoint(1);
%             StrainWithinEachElementGausspoint(2,4) = StrainEyyWithinEachElementNodalpoint(2);
%             StrainWithinEachElementGausspoint(3,4) = StrainEyyWithinEachElementNodalpoint(3);
%             StrainWithinEachElementGausspoint(4,4) = StrainEyyWithinEachElementNodalpoint(4);
% 
% 
%             % ------ Find the element nodal indices for strain ------
%             tempStrainIndex = [4*elements(j,1)-3  4*elements(j,1)-2  4*elements(j,1)-1  4*elements(j,1)  ...
%                                4*elements(j,2)-3  4*elements(j,2)-2  4*elements(j,2)-1  4*elements(j,2) ...
%                                4*elements(j,3)-3  4*elements(j,3)-2  4*elements(j,3)-1  4*elements(j,3) ...
%                                4*elements(j,4)-3  4*elements(j,4)-2  4*elements(j,4)-1  4*elements(j,4)];
% 
%             FStrain(tempStrainIndex) = FStrain(tempStrainIndex) + ...
%                [StrainWithinEachElementGausspoint(1,1); StrainWithinEachElementGausspoint(1,2); 
%                 StrainWithinEachElementGausspoint(1,3); StrainWithinEachElementGausspoint(1,4); 
%                 StrainWithinEachElementGausspoint(2,1); StrainWithinEachElementGausspoint(2,2); 
%                 StrainWithinEachElementGausspoint(2,3); StrainWithinEachElementGausspoint(2,4); 
%                 StrainWithinEachElementGausspoint(3,1); StrainWithinEachElementGausspoint(3,2); 
%                 StrainWithinEachElementGausspoint(3,3); StrainWithinEachElementGausspoint(3,4); 
%                 StrainWithinEachElementGausspoint(4,1); StrainWithinEachElementGausspoint(4,2); 
%                 StrainWithinEachElementGausspoint(4,3); StrainWithinEachElementGausspoint(4,4)];
% 
%             FStrainAvgTimes(tempStrainIndex) = FStrainAvgTimes(tempStrainIndex) + ones(16,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%% End of Comment %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
if waitBarDisplayOrNot == 0
    close(hbar);
end

% switch GaussPtOrder
%     case 2
%         Strain = FStrain./FStrainAvgTimes;
%     case 3
%         Strain = zeros(size(coordinates,1)*4,1);
% end

% ====== Use gridfit and Gauss points to fit strain field ======
Coordxnodes = [min(coordinatesFEM(:,1)):h:max(coordinatesFEM(:,1))]'; 
Coordynodes = [min(coordinatesFEM(:,2)):h:max(coordinatesFEM(:,2))]';
disp('** This may take a while when quadtree minimum element size is small. Please wait. ')
Iblur_10 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,1),Coordxnodes,Coordynodes); Iblur_10=Iblur_10';
disp('**** Please keep waiting! ')
Iblur_20 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,2),Coordxnodes,Coordynodes); Iblur_20=Iblur_20';
disp('****** Half way! Do not give up! ')
Iblur_30 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,3),Coordxnodes,Coordynodes); Iblur_30=Iblur_30';
disp('******** Almost there! ')
Iblur_40 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,4),Coordxnodes,Coordynodes); Iblur_40=Iblur_40';    
disp('********** Hoooooray!!! ********** ')
Strain = 0*[U(1:end-DIM*NodesNumPerEle); U(1:end-DIM*NodesNumPerEle)]; 

for tempi = 1:size(coordinatesFEM,1)
      [row1,~] = find(Coordxnodes == coordinatesFEM(tempi,1));
        [row2,~] = find(Coordynodes == coordinatesFEM(tempi,2));
        Strain(4*tempi-3) = Iblur_10(row1,row2);
        Strain(4*tempi-2) = Iblur_20(row1,row2);
        Strain(4*tempi-1) = Iblur_30(row1,row2);
        Strain(4*tempi-0) = Iblur_40(row1,row2);
end
    

end



%% ========= subroutines for  FEM Q4 shape function derivatives ========
function DN5Dksi=funDN5Dksi(ksi,eta,deltaPt5)
DN5Dksi = deltaPt5*0.5*(1-abs(eta)) ;
end
function DN5Deta=funDN5Deta(ksi,eta,deltaPt5)
DN5Deta = deltaPt5*0.5*(1+ksi)*sign(-eta);
end
function DN6Dksi=funDN6Dksi(ksi,eta,deltaPt6)
DN6Dksi = deltaPt6*0.5*(1+eta)*sign(-ksi);
end
function DN6Deta=funDN6Deta(ksi,eta,deltaPt6)
DN6Deta = deltaPt6*0.5*(1-abs(ksi));
end
function DN7Dksi=funDN7Dksi(ksi,eta,deltaPt7)
DN7Dksi = deltaPt7*0.5*(-1)*(1-abs(eta));
end
function DN7Deta=funDN7Deta(ksi,eta,deltaPt7)
DN7Deta = deltaPt7*0.5*(1-ksi)*sign(-eta);
end
function DN8Dksi=funDN8Dksi(ksi,eta,deltaPt8)
DN8Dksi = deltaPt8*0.5*(1-eta)*sign(-ksi);
end
function DN8Deta=funDN8Deta(ksi,eta,deltaPt8)
DN8Deta = deltaPt8*0.5*(-1)*(1-abs(ksi));
end
function DN1Dksi = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)
DN1Dksi = -0.25*(1-eta)-0.5*((deltaPt7*0.5*(-1)*(1-abs(eta)))+deltaPt8*0.5*(1-eta)*sign(-ksi));
end
function DN1Deta = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)
DN1Deta = -0.25*(1-ksi)-0.5*(deltaPt7*0.5*(1-ksi)*sign(-eta)+deltaPt8*0.5*(-1)*(1-abs(ksi)));
end
function DN2Dksi = funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)
DN2Dksi = 0.25*(1-eta)-0.5*(deltaPt8*0.5*(1-eta)*sign(-ksi)+deltaPt5*0.5*(1-abs(eta)));
end
function DN2Deta = funDN2Deta(ksi,eta,deltaPt8,deltaPt5)
DN2Deta = -0.25*(1+ksi)-0.5*(deltaPt8*0.5*(-1)*(1-abs(ksi))+deltaPt5*0.5*(1+ksi)*sign(-eta));
end
function DN3Dksi = funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)
DN3Dksi = 0.25*(1+eta)-0.5*(deltaPt5*0.5*(1-abs(eta))+deltaPt6*0.5*(1+eta)*sign(-ksi));
end
function DN3Deta = funDN3Deta(ksi,eta,deltaPt5,deltaPt6)
DN3Deta = 0.25*(1+ksi)-0.5*(deltaPt5*0.5*(1+ksi)*sign(-eta)+deltaPt6*0.5*(1-abs(ksi)));
end
function DN4Dksi = funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)
DN4Dksi = -0.25*(1+eta)-0.5*(deltaPt6*0.5*(1+eta)*sign(-ksi)+deltaPt7*0.5*(-1)*(1-abs(eta)));
end
function DN4Deta = funDN4Deta(ksi,eta,deltaPt6,deltaPt7)
DN4Deta = 0.25*(1-ksi)-0.5*(deltaPt6*0.5*(1-abs(ksi))+deltaPt7*0.5*(1-ksi)*sign(-eta));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 