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

contactTol = 5e-4; % 0.5 mm: CAD vertices closer than this across parts flag the shared edge/face
pad        = 1e-4; % 0.1 mm halo padded onto the localized zone box

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
% Locate the actual contact patch using each part's own CAD (B-rep)
% vertices -- no meshing, no scaling/growing of either part. Any vertex
% of A within contactTol of a vertex of B (and vice versa) sits on the
% shared edge/face loop; bound just those points, pad by a small halo,
% and subtract that small local box from BOTH parts. Nothing outside
% this box is touched.
    vA = A.Vertices;
    vB = B.Vertices;

    dA = minDistEach(vA, vB);
    closeA = vA(dA < contactTol, :);
    dB = minDistEach(vB, vA);
    closeB = vB(dB < contactTol, :);

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
    d = zeros(size(P,1), 1);
    for i = 1:size(P,1)
        diffs = Q - P(i,:);
        d(i) = sqrt(min(sum(diffs.^2, 2)));
    end
end
