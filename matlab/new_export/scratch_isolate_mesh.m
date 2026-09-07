clear; clc;
load('placed_parts.mat', 'placed', 'labels');

order = {'top_plate_8','bottom_plate_6','center_solenoid_3', ...
         'outer_solenoid_1','outer_solenoid_2','outer_solenoid_7','outer_solenoid_9', ...
         'chamber_5','injector_4'};

blob = placed{strcmp(labels, order{1})};
fprintf('%-16s alone: Faces=%d\n', order{1}, blob.NumFaces);
try
    m = generateMesh(blob, Hmax=0.01, Hmin=0.0005);
    fprintf('  [PASS] meshed: %d elements\n', size(m.Mesh.Elements,2));
catch ME
    fprintf('  [FAIL] %s\n', ME.message);
end

for i = 2:numel(order)
    part = placed{strcmp(labels, order{i})};
    blob = union(blob, part);
    fprintf('+ %-16s -> Cells=%d Faces=%d\n', order{i}, blob.NumCells, blob.NumFaces);
    try
        m = generateMesh(blob, Hmax=0.01, Hmin=0.0005);
        fprintf('  [PASS] meshed: %d elements\n', size(m.Mesh.Elements,2));
    catch ME
        fprintf('  [FAIL] %s\n', ME.message);
    end
end
