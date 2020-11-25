% Section 3* for quadtree mesh

%%
% ====== Remove finite elements where there is a hole ======

coordinatesFEMQuadtree = DICmesh.coordinatesFEM;
elementsFEMQuadtree = DICmesh.elementsFEM;
irregular = zeros(0,3);

% Define circle
C = ParHole(1:2); % circle center
R = ParHole(3); % circle radius
h = 2; % minimize element size in the refined quadtree mesh

% Generate a Quadtree mesh
while 1
    [~,mark4] = markCircle(coordinatesFEMQuadtree,[],elementsFEMQuadtree,C,R,h*2); % Don't delete this "*2"
    mark4 = find(mark4);
    [coordinatesFEMQuadtree,elementsFEMQuadtree,irregular] = QrefineR(coordinatesFEMQuadtree,elementsFEMQuadtree,irregular,mark4);
    if isempty(mark4)
        break
    end
end

% Remove elements within the center hole
[~,markOutside4] = markCircleInside(coordinatesFEMQuadtree,elementsFEMQuadtree,C,R);
elementsFEMQuadtree = elementsFEMQuadtree(markOutside4,:);
[markEleHoleEdge4,markEleFarOutside4] = markCircleInside(coordinatesFEMQuadtree,elementsFEMQuadtree,C,R+10);
markCoordHoleEdge = unique(elementsFEMQuadtree(markEleHoleEdge4,:));
if markCoordHoleEdge(1)==0, markCoordHoleEdge=markCoordHoleEdge(2:end); end

% Plot
clf; patch('Faces', elementsFEMQuadtree, 'Vertices', coordinatesFEMQuadtree, 'Facecolor','none','linewidth',1)
axis equal; axis tight; set(gca,'fontsize',20);

% Update mesh for considering hanging nodes
for tempj=1:size(irregular,1)
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,1:2), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,8)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,2:3), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,5)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,3:4), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,6)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[4,1]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,7)=irregular(tempj,3); end
    
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[2,1]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,8)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[3,2]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,5)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[4,3]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,6)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[1,4]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,7)=irregular(tempj,3); end
end


% Initialize variable U for the generated quadtree mesh
F_dispu = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),U0(1:2:end) );
F_dispv = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),U0(2:2:end) );

U0 = 0*coordinatesFEMQuadtree(:);
temp = F_dispu(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); U0(1:2:end)=temp(:);
temp = F_dispv(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); U0(2:2:end)=temp(:);

Plotdisp_show(U0,coordinatesFEMQuadtree,elementsFEMQuadtree(:,1:4));

DICmesh.coordinatesFEM = coordinatesFEMQuadtree;
DICmesh.elementsFEM = elementsFEMQuadtree;
DICmesh.elementMinSize = h;
DICmesh.coordinatesFEMWorld = [DICmesh.coordinatesFEM(:,1),size(ImgRefMask,1)+1-DICmesh.coordinatesFEM(:,2)];

