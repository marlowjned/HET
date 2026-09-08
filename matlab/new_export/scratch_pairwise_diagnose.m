%% scratch_pairwise_diagnose.m
% Diagnostic-only: for every pair of parts that actually touches (per the
% addCell-probing contact graph, same method as check_placed.m), union
% just that pair in isolation and check (a) whether the resulting face
% count is exactly additive (clean shared boundary) or inflated (sliver
% faces) and (b) whether the pair meshes at all. Isolating pairs this way
% localizes which specific mating interfaces are broken, instead of only
% learning "the accumulated blob fails somewhere." Does not modify or
% save any geometry.
clear; clc;
load('placed_parts.mat', 'placed', 'labels');
n = numel(placed);

allPts = cell2mat(cellfun(@(g) g.Vertices, placed, 'UniformOutput', false));
center = mean(allPts,1);
R = 8*max(vecnorm(allPts-center,2,2));

touching = {};
for i=1:n
    for j=i+1:n
        gm = fegeometry(multisphere(R));
        gm = translate(gm, center);
        ok1 = true; ok2 = true;
        try, gm = addCell(gm, placed{i}); catch, ok1=false; end
        if ok1
            try, gm = addCell(gm, placed{j}); catch, ok2=false; end
        end
        if ~ok1 || ~ok2
            touching(end+1,:) = {i,j}; %#ok<AGROW>
        end
    end
end
fprintf('Found %d touching pairs out of %d total.\n\n', size(touching,1), nchoosek(n,2));

fprintf('%-20s %-20s %-8s %-8s %-8s %-8s %s\n', ...
    'partA','partB','facesA','facesB','sum','union','mesh');
for k = 1:size(touching,1)
    i = touching{k,1}; j = touching{k,2};
    a = placed{i}; b = placed{j};
    sumAlone = a.NumFaces + b.NumFaces;
    try
        u = union(a,b);
        meshStr = 'n/a';
        try
            m = generateMesh(u, Hmax=0.01, Hmin=0.0005);
            meshStr = sprintf('PASS(%d el)', size(m.Mesh.Elements,2));
        catch ME2
            meshStr = ['FAIL: ' ME2.message];
        end
        fprintf('%-20s %-20s %-8d %-8d %-8d %-8d %s\n', ...
            labels{i}, labels{j}, a.NumFaces, b.NumFaces, sumAlone, u.NumFaces, meshStr);
    catch ME
        fprintf('%-20s %-20s %-8d %-8d %-8d UNION FAILED: %s\n', ...
            labels{i}, labels{j}, a.NumFaces, b.NumFaces, sumAlone, ME.message);
    end
end
