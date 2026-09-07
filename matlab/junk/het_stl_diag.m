%% het_stl_diag.m  — temporary STL import diagnostic, safe to delete after
clear; close all; clc;

stlDir = fullfile(fileparts(mfilename('fullpath')), 'stl');
files = dir(fullfile(stlDir, '*.stl'));

fprintf('Found %d STL files\n\n', numel(files));

%% Test 1: can we import each STL individually?
gms = cell(numel(files), 1);
for i = 1:numel(files)
    fpath = fullfile(stlDir, files(i).name);
    try
        gms{i} = fegeometry(fpath);
        fprintf('[PASS] %d: %s  (Cells=%d Faces=%d)\n', i, files(i).name, gms{i}.NumCells, gms{i}.NumFaces);
    catch ME
        fprintf('[FAIL] %d: %s — %s\n', i, files(i).name, ME.message);
        gms{i} = [];
    end
end

%% Test 2: add each STL part into an analytic air sphere
% addCell requires BOTH geometries to be fegeometry objects (multisphere
% returns a bare DiscreteGeometry, so it must be wrapped first), and each
% added cell must sit strictly inside one existing cell -- parts that touch
% or overlap a previously-added part will fail one at a time here, which is
% reported per-part instead of aborting so all conflicts surface at once.
fprintf('\n--- Combining into air sphere ---\n');
valid = find(~cellfun(@isempty, gms));
if isempty(valid)
    fprintf('[FAIL] No valid geometries to combine.\n');
else
    % Bounding radius/center from ALL parts' vertices, not just the first
    allPts = cell2mat(cellfun(@(g) g.Vertices, gms(valid), 'UniformOutput', false));
    center = mean(allPts, 1);
    roughRadius = 6 * max(vecnorm(allPts - center, 2, 2));
    gm = fegeometry(multisphere(roughRadius));
    gm = translate(gm, center);

    added = false(numel(gms), 1);
    for i = valid'
        try
            gm = addCell(gm, gms{i});
            added(i) = true;
            fprintf('[PASS] added %s -> NumCells=%d\n', files(i).name, gm.NumCells);
        catch ME
            fprintf('[FAIL] could not add %s — %s\n', files(i).name, ME.message);
        end
    end

    fprintf('\nCombined geometry: Cells=%d Faces=%d (air + %d/%d parts added)\n', ...
        gm.NumCells, gm.NumFaces, nnz(added), numel(valid));
    if nnz(added) < numel(valid)
        fprintf(['Parts that failed to add likely touch or overlap an already-added part\n' ...
            '(shared/coincident faces from CAD mating). addCell needs strict interior\n' ...
            'placement -- union same-material touching bodies in Onshape before export,\n' ...
            'or add a small clearance gap, then re-export and re-run.\n']);
    end

    figure('Name', 'STL combined in air sphere');
    pdegplot(gm, CellLabels='on', FaceAlpha=0.1);
    title('Parts in air sphere — rotate to identify each cell');
end
