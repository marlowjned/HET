clear; clc;
load('placed_parts.mat', 'placed', 'labels');

R = [1 0 0; 0 0 -1; 0 1 0];  % from T table: Zax=[0 -1 0], Xax=[1 0 0], Yax=cross=[0 0 1]... recompute properly below
Zax = [0 -1 0]; Xax = [1 0 0]; Yax = cross(Zax,Xax);
R = [Xax(:) Yax(:) Zax(:)];
loc = [-0.0615335069596768, 0.0798566981732845, 0.0135308774188161];

chamber = placed{strcmp(labels,'chamber_5')};
% Transform chamber mesh nodes into the center_solenoid pole's local frame
gmC = generateMesh(chamber, Hmax=0.004);
p = gmC.Mesh.Nodes';                 % Nx3 global
pLocal = (R' * (p - loc)')';         % Nx3 local
rLocal = sqrt(pLocal(:,1).^2 + pLocal(:,2).^2);
zLocal = pLocal(:,3);

% Bin by z and report min radius (inner bore) of chamber material in each bin,
% restricted to the pole's own axial span [0, 0.11113]
edges = linspace(0, 0.11113, 12);
fprintf('Chamber radial extent vs local z (pole axis frame), pole radius=0.0308:\n');
for i=1:numel(edges)-1
    inBin = zLocal >= edges(i) & zLocal < edges(i+1);
    if any(inBin)
        fprintf('  z=[%.4f %.4f]: chamber r range=[%.4f %.4f]  (n=%d)\n', ...
            edges(i), edges(i+1), min(rLocal(inBin)), max(rLocal(inBin)), nnz(inBin));
    else
        fprintf('  z=[%.4f %.4f]: no chamber material here\n', edges(i), edges(i+1));
    end
end
