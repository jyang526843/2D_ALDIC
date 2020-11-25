% ==================================
% function compute strain Quadtree
% ----------------------------------
% switch MethodToComputeStrain
%   case 3: finite element Gauss points;
% ==================================
 

%% Update infinitesimal strain to large deformation gradient tensor
FStrain = FSubpb2;
FStrainFinite = FStrain;
for tempi = 1:4:length(FStrain)
    
    % Obtain each component of def grad tensor
    dudx = FStrain(tempi);
    dvdx = FStrain(tempi+1);
    dudy = FStrain(tempi+2);
    dvdy = FStrain(tempi+3);
    
    switch DICpara.StrainType
        case 0 % Infinitesimal stran
            % Do nothing
        case 1 % Eluerian strain
            FStrainFinite(tempi) = 1/(1-dudx)-1;
            FStrainFinite(tempi+3) = 1/(1-dvdy)-1;
            FStrainFinite(tempi+2) = dudy/(1-dvdy);
            FStrainFinite(tempi+1) = dvdx/(1-dudx);
        case 2 % Green-Lagrangian strain: E=(C-I)/2
            FStrainFinite(tempi) = 0.5*(dudx*2-dudx^2-dvdx^2);
            FStrainFinite(tempi+3) = 0.5*(dvdy*2-dudy^2-dvdy^2);
            FStrainFinite(tempi+2) = 0.5*(dudy+dvdx-dudx*dudy-dvdx*dvdy);
            FStrainFinite(tempi+1) = 0.5*(dvdx+dudy-dudy*dudx-dvdy*dvdx);
        case 3
            disp('Press "Ctrl+C" to modify by yourself.'); pause;
        otherwise
            disp('Wrong strain type!');
    end

end

FStraintemp = FStrainFinite;



