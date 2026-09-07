load('placed_parts.mat', 'placed', 'labels');
ironIdx = find(contains(labels, {'top_plate','bottom_plate','center_solenoid','outer_solenoid'}));
ironBlob = placed{ironIdx(1)};
for k = 2:numel(ironIdx)
    ironBlob = union(ironBlob, placed{ironIdx(k)});
end
chamberSolid = placed{strcmp(labels,'chamber_5')};

modes = {'default', 'keepboth'};
for m_ = 1:2
    if strcmp(modes{m_}, 'default')
        blob = union(ironBlob, chamberSolid);
    else
        blob = union(ironBlob, chamberSolid, KeepBoundaries=[true true]);
    end
    fprintf('%-10s Cells=%d Faces=%d\n', modes{m_}, blob.NumCells, blob.NumFaces);
    try
        mesh = generateMesh(blob, Hmax=0.008, Hmin=0.0004);
        fprintf('  [PASS] meshed: %d elements\n', size(mesh.Mesh.Elements,2));
    catch ME
        fprintf('  [FAIL] %s\n', ME.message);
    end
end
