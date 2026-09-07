ironFiles = {'BPL-700 Assembly - top_plate-1__Body1.step', ...
             'BPL-700 Assembly - bottom_plate-1__Body1.step', ...
             'BPL-700 Assembly - center_solenoid-1__Body1.step', ...
             'BPL-700 Assembly - outer_solenoid-1__Body1.step'};
iron = fegeometry(ironFiles{1});
for i=2:numel(ironFiles)
    iron = union(iron, fegeometry(ironFiles{i}));
end

chamber  = fegeometry('BPL-700 Assembly - chamber-1__Body1.step');
injector = fegeometry('BPL-700 Assembly - injector-1__Body1.step');

contactTol = 2e-4; % 0.2 mm: nodes closer than this across parts are "at the interface"
pad        = 5e-5; % 50 micron halo padded onto the localized zone box

[ironG, chamberG]     = gapAtInterface(iron,    chamber,  contactTol, pad);
[ironG, injectorG]    = gapAtInterface(ironG,   injector, contactTol, pad);
[chamberG, injectorG] = gapAtInterface(chamberG,injectorG,contactTol, pad);

fprintf('iron    NumCells after gapping (should stay 1): %d\n', ironG.NumCells);
fprintf('chamber NumCells after gapping (should stay 1): %d\n', chamberG.NumCells);
fprintf('injector NumCells after gapping (should stay 1): %d\n', injectorG.NumCells);

allPts = [ironG.Vertices; chamberG.Vertices; injectorG.Vertices];
center = mean(allPts,1);
R = 8*max(vecnorm(allPts-center,2,2));
gm = fegeometry(multisphere(R));
gm = translate(gm, center);
parts = {ironG, chamberG, injectorG};
labels = {'iron','chamber','injector'};
for i=1:numel(parts)
    try
        gm = addCell(gm, parts{i});
        fprintf('[PASS] added %s -> NumCells=%d\n', labels{i}, gm.NumCells);
    catch ME
        fprintf('[FAIL] %s -- %s\n', labels{i}, ME.message);
    end
end

function [A2,B2] = gapAtInterface(A,B,contactTol,pad)
% Locate the actual contact patch by meshing both parts and finding mesh
% nodes on A within contactTol of any node on B (and vice versa) -- this
% is the true interface locus, not a whole-part bounding box. Bound just
% those near-contact nodes, pad by a small halo, and subtract that small
% local box from BOTH parts. Nothing outside this box is touched.
    Am = generateMesh(A, Hmax=contactTol*15);
    Bm = generateMesh(B, Hmax=contactTol*15);
    nodesA = Am.Mesh.Nodes';   % Nx3
    nodesB = Bm.Mesh.Nodes';

    dA = minDistEach(nodesA, nodesB);
    closeA = nodesA(dA < contactTol, :);
    dB = minDistEach(nodesB, nodesA);
    closeB = nodesB(dB < contactTol, :);

    contactPts = [closeA; closeB];
    if isempty(contactPts)
        A2 = A; B2 = B;   % not actually touching within tolerance
        return;
    end

    lo = min(contactPts,[],1) - pad;
    hi = max(contactPts,[],1) + pad;
    dims = hi - lo;
    ctr  = (lo + hi) / 2;
    zone = fegeometry(multicuboid(dims(1), dims(2), dims(3)));
    zone = translate(zone, ctr - [0 0 dims(3)/2]);
    A2 = subtract(A, zone);
    B2 = subtract(B, zone);
end

function d = minDistEach(P, Q)
% No Statistics Toolbox available -- plain vectorized nearest-neighbor
% distance from each row of P to the closest row of Q.
    d = zeros(size(P,1), 1);
    for i = 1:size(P,1)
        diffs = Q - P(i,:);
        d(i) = sqrt(min(sum(diffs.^2, 2)));
    end
end
