%% het_view_assembly.m
%
% Quick visual sanity check of the reassembled BPL-700 geometry: loads the
% correctly-placed parts (positions/orientations recovered from the real
% Onshape assembly transforms -- see new_export/apply_transforms.m) and
% plots them together, color-coded by part, with centroid labels so you
% can rotate the figure and confirm nothing is mispositioned, duplicated,
% or still colliding before trusting it for the B-field solve.
%
% Requires new_export/placed_parts.mat (run apply_transforms.m first).

clear; close all; clc;

matFile = fullfile(fileparts(mfilename('fullpath')), 'new_export', 'placed_parts.mat');
load(matFile, 'placed', 'labels');

% One color per physical part identity (not per occurrence), so the 4
% outer_solenoid copies read as "the same part in 4 places."
baseName = regexprep(labels, '_\d+$', '');
uniqueNames = unique(baseName, 'stable');
palette = lines(numel(uniqueNames));
colorOf = containers.Map(uniqueNames, num2cell(palette, 2));

figure('Name', 'BPL-700 reassembled geometry', 'Color', 'w');

% pdegplot resets "hold" to off internally after every call, so it must be
% re-asserted before each individual call or later parts wipe out earlier
% ones instead of adding to the plot.
for i = 1:numel(placed)
    nm = baseName{i};
    c = colorOf(nm);

    hold on;
    pdegplot(placed{i}, FaceAlpha=0.55, FaceColor=c, Lighting='off');

    ctr = mean(placed{i}.Vertices, 1);
    hold on;
    text(ctr(1), ctr(2), ctr(3), sprintf('  %s', labels{i}), ...
        FontSize=8, Color=[0 0 0], FontWeight='bold');
end

% Dummy proxy markers purely for a reliable legend, added after everything
% else so they aren't wiped by a later pdegplot/hold reset.
hold on;
legendHandles = gobjects(numel(uniqueNames), 1);
for k = 1:numel(uniqueNames)
    legendHandles(k) = plot3(NaN, NaN, NaN, 's', ...
        MarkerFaceColor=colorOf(uniqueNames{k}), MarkerEdgeColor='none', MarkerSize=10);
end

axis equal;
grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('BPL-700 reassembled geometry (rotate to inspect)');
legend(legendHandles, uniqueNames, Location='eastoutside', Interpreter='none');
view(35, 20);
camlight('headlight');
lighting gouraud;

fprintf('Plotted %d parts (%d unique identities).\n', numel(placed), numel(uniqueNames));
fprintf('Rotate/zoom the figure to check: correct part count, no stray duplicates,\n');
fprintf('plausible touching, and sane overall proportions.\n');
