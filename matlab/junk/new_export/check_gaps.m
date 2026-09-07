load('placed_parts.mat', 'placed', 'labels');
bpIdx = find(strcmp(labels,'bottom_plate_6'));
tpIdx = find(strcmp(labels,'top_plate_8'));
bp = placed{bpIdx};
tp = placed{tpIdx};
for i=1:numel(labels)
    if startsWith(labels{i}, 'outer_solenoid')
        os = placed{i};
        dBP = min(min(pdist2local(os.Vertices, bp.Vertices)));
        dTP = min(min(pdist2local(os.Vertices, tp.Vertices)));
        fprintf('%-20s min-dist to bottom_plate = %.6g m,  to top_plate = %.6g m\n', labels{i}, dBP, dTP);
    end
end

function D = pdist2local(A,B)
    D = zeros(size(A,1), size(B,1));
    for i=1:size(A,1)
        D(i,:) = sqrt(sum((B - A(i,:)).^2, 2))';
    end
end
