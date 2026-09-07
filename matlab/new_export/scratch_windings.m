clear; clc;
load('placed_parts.mat', 'placed', 'labels');

% Same transform table as apply_transforms.m, but we only need the 5 pole
% occurrences here (center_solenoid + outer_solenoid x4).
T = {
  'outer_solenoid',  [-0.125033506959677, -0.0217433018267154,  0.0770308774188161], [0 1 0], [1 0 0]
  'outer_solenoid',  [ 0.00196649304032328, 0.0798566981732845,  0.0770308774188161], [0 -1 0], [1 0 0]
  'center_solenoid', [-0.0615335069596768,  0.0798566981732845,  0.0135308774188161], [0 -1 0], [1 0 0]
  'outer_solenoid',  [ 0.00196649304032332, 0.0798566981732845, -0.0499691225811839], [0 -1 0], [1 0 0]
  'outer_solenoid',  [-0.125033506959677,  -0.0217433018267154, -0.0499691225811839], [0 1 0],  [1 0 0]
};
poleLabels = {'outer_solenoid_1','outer_solenoid_2','center_solenoid_3','outer_solenoid_7','outer_solenoid_9'};

% Winding dimensions per pole type, in LOCAL frame -- explicit id/od/z0/z1,
% not a generic formula, because center_solenoid is a stepped pole (thin
% shaft for most of its length, flush/thick only in the last ~11mm against
% the chamber -- see scratch_center_gap.m / profiling). The winding must
% sit on the THIN section, using the chamber's own bore (~30.6mm radius)
% as the natural OD limit; the outer poles are uniform cylinders with
% open radial clearance everywhere, so a simple clearance+wall works.
clearance = 0.002; % 2 mm radial clearance from any neighboring solid

outer.id = 2*(0.01588 + clearance);
outer.od = outer.id + 2*0.008;         % 8 mm wall, arbitrary but reasonable
outer.z0 = 0.00762;                    % inset from both pole ends (~15% margin)
outer.z1 = outer.z0 + 0.10160*0.85;

center.id = 2*(0.01588 + clearance);   % thin-section pole radius + clearance
center.od = 2*(0.0306 - clearance);    % chamber bore radius - clearance
center.z0 = 0.005;                     % margin from bottom_plate end (z=0)
center.z1 = 0.096;                     % margin before the step to the thick section (~0.101)

dims.outer_solenoid  = outer;
dims.center_solenoid = center;

windings = cell(size(T,1),1);
Rlist    = cell(size(T,1),1);
loclist  = cell(size(T,1),1);

for i = 1:size(T,1)
    name = T{i,1};
    loc  = T{i,2};
    Zax  = T{i,3} / norm(T{i,3});
    Xax  = T{i,4} / norm(T{i,4});
    Yax  = cross(Zax, Xax);
    R = [Xax(:), Yax(:), Zax(:)];

    d = dims.(name);
    h = d.z1 - d.z0;

    sleeve = fegeometry(multicylinder([d.id/2, d.od/2], h, Void=[true false]));
    sleeve = translate(sleeve, [0 0 d.z0]);   % local frame placement along pole axis

    [axisVec, angleDeg] = rotmat2axisangle(R);
    if angleDeg > 1e-6
        sleeve = rotate(sleeve, angleDeg, [0 0 0], axisVec);
    end
    sleeve = translate(sleeve, loc);

    windings{i} = sleeve;
    Rlist{i} = R;
    loclist{i} = loc;
    fprintf('winding for %-16s centroid=[%8.4f %8.4f %8.4f]  ID=%.4f OD=%.4f h=%.4f\n', ...
        poleLabels{i}, mean(sleeve.Vertices,1), d.id, d.od, h);
end

save('windings.mat', 'windings', 'poleLabels', 'Rlist', 'loclist');

%% Verify: windings should float freely -- not touch iron, chamber, injector, or each other
allPts = cell2mat(cellfun(@(g) g.Vertices, placed, 'UniformOutput', false));
center_ = mean(allPts,1);
Rair = 8*max(vecnorm(allPts-center_,2,2));

ironIdx = find(contains(labels, {'top_plate','bottom_plate','center_solenoid','outer_solenoid'}));
chamberIdx = find(strcmp(labels,'chamber_5'));
injectorIdx = find(strcmp(labels,'injector_4'));

ironBlob = placed{ironIdx(1)};
for k = 2:numel(ironIdx)
    ironBlob = union(ironBlob, placed{ironIdx(k)});
end
structureBlob = union(union(ironBlob, placed{chamberIdx}), placed{injectorIdx});
fprintf('\nstructureBlob NumCells (should be 1): %d\n', structureBlob.NumCells);

gm = fegeometry(multisphere(Rair));
gm = translate(gm, center_);
gm = addCell(gm, structureBlob);
fprintf('after adding structureBlob: NumCells=%d\n', gm.NumCells);

for i = 1:numel(windings)
    try
        gm = addCell(gm, windings{i});
        fprintf('[PASS] added winding %d (%s) -> NumCells=%d\n', i, poleLabels{i}, gm.NumCells);
    catch ME
        fprintf('[FAIL] winding %d (%s) -- %s\n', i, poleLabels{i}, ME.message);
    end
end

% Also verify windings don't collide with each other
for i = 1:numel(windings)
    for j = i+1:numel(windings)
        gm2 = fegeometry(multisphere(Rair));
        gm2 = translate(gm2, center_);
        ok = true;
        try
            gm2 = addCell(gm2, windings{i});
            gm2 = addCell(gm2, windings{j});
        catch
            ok = false;
        end
        if ~ok
            fprintf('[FAIL] winding %d and %d overlap each other\n', i, j);
        end
    end
end

function [axisVec, angleDeg] = rotmat2axisangle(R)
    c = max(-1, min(1, (trace(R)-1)/2));
    angleRad = acos(c);
    if angleRad < 1e-9
        axisVec = [1 0 0]; angleDeg = 0; return;
    end
    if abs(pi - angleRad) < 1e-6
        M = (R + eye(3)) / 2;
        [~, idx] = max(diag(M));
        v = sqrt(max(M(:,idx), 0))';
        others = setdiff(1:3, idx);
        for k = others
            if M(idx,k) < 0
                v(k) = -v(k);
            end
        end
        axisVec = v / norm(v);
        angleDeg = 180;
        return;
    end
    v = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)] / (2*sin(angleRad));
    axisVec = (v / norm(v))';
    angleDeg = rad2deg(angleRad);
end
