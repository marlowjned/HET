%% scratch_check_tessellation.m
% Diagnostic-only: checks each part's own STEP-derived boundary
% triangulation (independent of any union) for degenerate near-zero-area
% triangles. This is a per-part export-quality check, distinct from the
% mating-interface check in scratch_pairwise_diagnose.m -- a part can be
% individually clean and still produce a bad union with its neighbor (a
% placement/tolerance problem), or individually dirty regardless of who
% it's unioned with (an export-tessellation problem). Does not modify or
% save any geometry.
clear; clc;
load('placed_parts.mat', 'placed', 'labels');

fprintf('%-20s %-8s %-14s %-14s %s\n', 'part','faces','minTriArea','#tri<1e-10','#tri<1e-8');
for i = 1:numel(labels)
    g = placed{i};
    try
        TR = triangulation(g);
        P = TR.Points; C = TR.ConnectivityList;
        v1 = P(C(:,2),:) - P(C(:,1),:);
        v2 = P(C(:,3),:) - P(C(:,1),:);
        areas = 0.5*vecnorm(cross(v1,v2,2),2,2);
        fprintf('%-20s %-8d %-14.3e %-14d %d\n', ...
            labels{i}, g.NumFaces, min(areas), nnz(areas<1e-10), nnz(areas<1e-8));
    catch ME
        fprintf('%-20s triangulation FAILED: %s\n', labels{i}, ME.message);
    end
end
