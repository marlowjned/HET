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

pad = 5e-5; % 50 micron halo around the bbox-overlap contact patch

[ironG, chamberG]     = gapAtInterface(iron,    chamber,  pad);
[ironG, injectorG]    = gapAtInterface(ironG,   injector, pad);
[chamberG, injectorG] = gapAtInterface(chamberG,injectorG,pad);

fprintf('iron   : orig verts=%d -> gapped verts=%d\n', size(iron.Vertices,1), size(ironG.Vertices,1));
fprintf('chamber: orig verts=%d -> gapped verts=%d\n', size(chamber.Vertices,1), size(chamberG.Vertices,1));
fprintf('injector: orig verts=%d -> gapped verts=%d\n', size(injector.Vertices,1), size(injectorG.Vertices,1));

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

function [A2,B2] = gapAtInterface(A,B,pad)
% Localizes the contact patch as the overlap of A's and B's bounding
% boxes (their true shared face must lie in this region), pads it by a
% small halo, and subtracts that small box from BOTH parts. Does not
% touch either part's geometry anywhere outside this local box.
    bbA = [min(A.Vertices,[],1); max(A.Vertices,[],1)];
    bbB = [min(B.Vertices,[],1); max(B.Vertices,[],1)];
    lo = max(bbA(1,:), bbB(1,:)) - pad;
    hi = min(bbA(2,:), bbB(2,:)) + pad;
    if any(lo >= hi)
        A2 = A; B2 = B;   % bounding boxes don't overlap -- not touching
        return;
    end
    dims = hi - lo;
    ctr  = (lo + hi) / 2;
    zone = fegeometry(multicuboid(dims(1), dims(2), dims(3)));
    zone = translate(zone, ctr - [0 0 dims(3)/2]);  % multicuboid: centered x/y, z in [0,dims(3)]
    A2 = subtract(A, zone);
    B2 = subtract(B, zone);
end
