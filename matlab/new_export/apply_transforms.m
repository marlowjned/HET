%% apply_transforms.m
% Applies the real per-occurrence placement transforms (extracted from the
% whole-assembly STEP export BPL-700-assem-SIMS-V1.step) to the reliable
% flat per-part STEP geometries, fixing the outer_solenoid duplicate-
% position bug and putting every part in its true assembled location.

clear; clc;

% name, translation (loc), Z-axis direction, X-refdir -- extracted from
% NAUO->IDT chain in BPL-700-assem-SIMS-V1.step (see scratch_check_transforms.py)
T = {
  'outer_solenoid',  [-0.125033506959677, -0.0217433018267154,  0.0770308774188161], [0 1 0], [1 0 0]
  'outer_solenoid',  [ 0.00196649304032328, 0.0798566981732845,  0.0770308774188161], [0 -1 0], [1 0 0]
  'center_solenoid', [-0.0615335069596768,  0.0798566981732845,  0.0135308774188161], [0 -1 0], [1 0 0]
  'injector',        [-0.0615926647537434,  0.0766816981732846,  0.0135553813794385], [0 -1 0], [1 0 0]
  'chamber',         [-0.0615926647537434,  0.0798566981732845,  0.0135553813794386], [0 -1 0], [1 0 0]
  'bottom_plate',    [-0.0615335069596767,  0.0893816981732846,  0.0135308774188161], [0 0 1],  [1 0 0]
  'outer_solenoid',  [ 0.00196649304032332, 0.0798566981732845, -0.0499691225811839], [0 -1 0], [1 0 0]
  'top_plate',       [-0.0615335069596767, -0.0217433018267154,  0.0135308774188161], [0 0 1],  [1 0 0]
  'outer_solenoid',  [-0.125033506959677,  -0.0217433018267154, -0.0499691225811839], [0 1 0],  [1 0 0]
};

fileFor = containers.Map( ...
    {'top_plate','bottom_plate','chamber','injector','center_solenoid','outer_solenoid'}, ...
    {'BPL-700 Assembly - top_plate-1__Body1.step', ...
     'BPL-700 Assembly - bottom_plate-1__Body1.step', ...
     'BPL-700 Assembly - chamber-1__Body1.step', ...
     'BPL-700 Assembly - injector-1__Body1.step', ...
     'BPL-700 Assembly - center_solenoid-1__Body1.step', ...
     'BPL-700 Assembly - outer_solenoid-1__Body1.step'});

placed = cell(size(T,1),1);
labels = cell(size(T,1),1);
for i = 1:size(T,1)
    name = T{i,1};
    loc  = T{i,2};
    Zax  = T{i,3} / norm(T{i,3});
    Xax  = T{i,4} / norm(T{i,4});
    Yax  = cross(Zax, Xax);
    R = [Xax(:), Yax(:), Zax(:)];   % columns = global directions of local X,Y,Z

    gm = fegeometry(fileFor(name));

    [axisVec, angleDeg] = rotmat2axisangle(R);
    if angleDeg > 1e-6
        gm = rotate(gm, angleDeg, [0 0 0], axisVec);
    end
    gm = translate(gm, loc);

    placed{i} = gm;
    labels{i} = sprintf('%s_%d', name, i);
    fprintf('%-20s centroid=[%8.4f %8.4f %8.4f]\n', labels{i}, mean(gm.Vertices,1));
end

save('placed_parts.mat', 'placed', 'labels');
fprintf('\nSaved placed_parts.mat with %d correctly-positioned parts.\n', numel(placed));

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
        % fix signs via off-diagonal terms
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
