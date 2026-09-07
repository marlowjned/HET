clear; clc;
load('placed_parts.mat', 'placed', 'labels');
n = numel(placed);

allPts = cell2mat(cellfun(@(g) g.Vertices, placed, 'UniformOutput', false));
center = mean(allPts,1);
fprintf('Overall assembly bbox: x=[%.4f %.4f] y=[%.4f %.4f] z=[%.4f %.4f]\n', ...
    min(allPts(:,1)), max(allPts(:,1)), min(allPts(:,2)), max(allPts(:,2)), ...
    min(allPts(:,3)), max(allPts(:,3)));

fprintf('\nPairwise touch/overlap check (%d pairs):\n', nchoosek(n,2));
R = 8*max(vecnorm(allPts-center,2,2));
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
            fprintf('  TOUCH: %-20s <-> %-20s\n', labels{i}, labels{j});
        end
    end
end
