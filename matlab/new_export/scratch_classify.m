clear; clc;
load('placed_parts.mat', 'placed', 'labels');
load('windings.mat', 'windings', 'poleLabels');

ironIdx     = find(contains(labels, {'top_plate','bottom_plate','center_solenoid','outer_solenoid'}));
chamberIdx  = find(strcmp(labels,'chamber_5'));
injectorIdx = find(strcmp(labels,'injector_4'));

ironBlob = placed{ironIdx(1)};
for k = 2:numel(ironIdx)
    ironBlob = union(ironBlob, placed{ironIdx(k)});
end
chamberSolid  = placed{chamberIdx};
injectorSolid = placed{injectorIdx};
structureBlob = union(union(ironBlob, chamberSolid), injectorSolid);

allPts = cell2mat(cellfun(@(g) g.Vertices, placed, 'UniformOutput', false));
center = mean(allPts,1);
airRadius = 7*max(vecnorm(allPts-center,2,2));

gm = fegeometry(multisphere(airRadius));
gm = translate(gm, center);
gm = addCell(gm, structureBlob);          % cell 2 = structureBlob
for i = 1:numel(windings)
    gm = addCell(gm, windings{i});        % cells 3..7 = windings
end
fprintf('Assembled geometry: NumCells=%d\n', gm.NumCells);

% Coarse-ish mesh for a first mechanical check of the classify pipeline.
gm = generateMesh(gm, Hmax=airRadius/10, Hface={cellFaces(gm,2:gm.NumCells), 0.006});
fprintf('Mesh: %d elements, %d nodes\n', size(gm.Mesh.Elements,2), size(gm.Mesh.Nodes,2));

nElem = size(gm.Mesh.Elements,2);
regionID = zeros(nElem,1);

% Air stays region 1
idx = findElements(gm.Mesh, 'region', Cell=1);
regionID(idx) = 1;

% Windings: cell k (3..7) -> region (k+2), i.e. 5..9
for k = 3:gm.NumCells
    idx = findElements(gm.Mesh, 'region', Cell=k);
    regionID(idx) = k + 2;
end

% structureBlob (cell 2): classify by which original solid each element
% centroid falls inside
idxBlob = findElements(gm.Mesh, 'region', Cell=2);
elemNodes = gm.Mesh.Elements(1:4, idxBlob);   % first 4 = corner nodes even for quadratic elements
coords = gm.Mesh.Nodes;
centroids = squeeze(mean(reshape(coords(:, elemNodes), 3, 4, []), 2))';  % Nx3

inIron     = ~isnan(findCell(ironBlob, centroids));
inChamber  = ~isnan(findCell(chamberSolid, centroids));
inInjector = ~isnan(findCell(injectorSolid, centroids));

fprintf('structureBlob elements: %d total | iron=%d chamber=%d injector=%d unclassified=%d\n', ...
    numel(idxBlob), nnz(inIron), nnz(inChamber), nnz(inInjector), ...
    numel(idxBlob) - nnz(inIron|inChamber|inInjector));

regionID(idxBlob(inIron))     = 2;
regionID(idxBlob(inChamber))  = 3;
regionID(idxBlob(inInjector)) = 4;

unresolved = idxBlob(~(inIron|inChamber|inInjector));
if ~isempty(unresolved)
    fprintf('Resolving %d unclassified elements by nearest solid...\n', numel(unresolved));
    for ii = 1:numel(unresolved)
        e = unresolved(ii);
        p = mean(coords(:, elemNodes(:, idxBlob==e))', 1);
        dIron = min(vecnorm(ironBlob.Vertices - p, 2, 2));
        dCh   = min(vecnorm(chamberSolid.Vertices - p, 2, 2));
        dInj  = min(vecnorm(injectorSolid.Vertices - p, 2, 2));
        [~, which] = min([dIron, dCh, dInj]);
        regionID(e) = which + 1;  % 2,3,4
    end
end

fprintf('\nFinal region element counts:\n');
for r = 1:9
    fprintf('  region %d: %d elements\n', r, nnz(regionID==r));
end

save('classified.mat', 'gm', 'regionID', '-v7.3');
