%% het_cad_inspect.m
%
% First-pass tool for identifying unlabeled solid bodies in an imported
% STEP assembly: plots cell labels and prints a per-cell geometric
% fingerprint (bounding box, centroid, size) so cell index -> physical
% part (core, pole piece, coil, housing, ...) can be matched by hand.
% Output feeds the cell-ID lookup table used by the CAD-based solve script.
%
% Requires MATLAB R2023b+ with PDE Toolbox.

clear; close all; clc;

stepFile = fullfile(fileparts(mfilename('fullpath')), 'BPL-700-part-chamber.step');
airDomainFactor = 6; % air sphere radius = this x max distance from assembly centroid

%% Import + recenter
cadGm = fegeometry(stepFile);
center = mean(cadGm.Vertices, 1);
cadGm = translate(cadGm, -center);

fprintf('Imported %s\n', stepFile);
fprintf('Cells: %d, Faces: %d, original centroid: [%.4f %.4f %.4f] m\n\n', ...
    cadGm.NumCells, cadGm.NumFaces, center);

%% Wrap in air domain
maxExtent = max(vecnorm(cadGm.Vertices - center, 2, 2));
airRadius = airDomainFactor * maxExtent;
gm = multisphere(airRadius);
gm = addCell(gm, cadGm);
cadCellRange = 2:gm.NumCells; % cell 1 is air; CAD cells appended after

%% Thrust-axis guess: for a squat disc-shaped thruster, the thrust axis
% usually has the smallest bounding-box extent (radial directions are
% wider). Heuristic only -- confirm against the labeled plot below.
ranges = max(cadGm.Vertices, [], 1) - min(cadGm.Vertices, [], 1);
axisNames = {'X', 'Y', 'Z'};
[~, likelyAxisIdx] = min(ranges);
fprintf('Likely thrust axis (smallest extent): %s  (X=%.4f Y=%.4f Z=%.4f m)\n\n', ...
    axisNames{likelyAxisIdx}, ranges(1), ranges(2), ranges(3));

figure('Name', 'CAD cell labels');
pdegplot(gm, CellLabels="on", FaceAlpha=0.1);
hold on;
axLen = 1.2 * maxExtent;
line([0 axLen], [0 0], [0 0], Color="r", LineWidth=2); text(axLen, 0, 0, 'X', Color="r");
line([0 0], [0 axLen], [0 0], Color="g", LineWidth=2); text(0, axLen, 0, 'Y', Color="g");
line([0 0], [0 0], [0 axLen], Color="b", LineWidth=2); text(0, 0, axLen, 'Z', Color="b");
hold off;
title('BPL-700 assembly: identify each cell number against the CAD; confirm thrust axis (colored lines)');

%% Coarse mesh, just to extract per-cell node positions (no solve)
model = femodel(AnalysisType="magnetostatic", Geometry=gm);
model.MaterialProperties = materialProperties(RelativePermeability=1);
model = generateMesh(model, Hmax=airRadius/8);

fprintf('%-6s %10s %10s %28s %28s\n', 'Cell', 'NumElem', 'Volume*', 'BBox size [x y z] (m)', 'Centroid [x y z] (m)');
for c = cadCellRange
    elemIdx = findElements(model.Mesh, "region", Cell=c);
    nodeIdx = unique(model.Mesh.Elements(:, elemIdx));
    pts = model.Mesh.Nodes(:, nodeIdx)';
    bboxSize = max(pts, [], 1) - min(pts, [], 1);
    centroid = mean(pts, 1);
    fprintf('%-6d %10d %10s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', ...
        c, numel(elemIdx), '-', bboxSize, centroid);
end
fprintf('\n*Volume omitted (node-based bbox only); use bbox size + centroid to match parts.\n');
fprintf('Rotate the figure and cross-reference cell numbers against the Onshape assembly.\n');
