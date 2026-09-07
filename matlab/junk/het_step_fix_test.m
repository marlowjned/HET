%% het_step_fix_test.m — test the Python STEP rewriter, safe to delete after
clear; clc;
base = fileparts(mfilename('fullpath'));

files = {
    'BPL-700-part-chamber.step'
    'BPL-700-sims-assem.step'
    'BPL-700-sims-assem-203.step'
};

for i = 1:numel(files)
    inFile  = fullfile(base, files{i});
    [~,n,~] = fileparts(inFile);
    outFile = fullfile(base, [n '_clean.step']);

    fprintf('=== %s ===\n', files{i});

    % Run the Python rewriter
    cmd = sprintf('python "%s" "%s"', fullfile(base,'fix_onshape_step.py'), inFile);
    [status, out] = system(cmd);
    fprintf('%s\n', out);
    if status ~= 0
        fprintf('[ERROR] Python script failed (exit %d)\n\n', status);
        continue;
    end

    % Try to import the cleaned file
    try
        gm = fegeometry(outFile);
        fprintf('[PASS] fegeometry OK — Cells=%d Faces=%d\n\n', gm.NumCells, gm.NumFaces);
    catch ME
        fprintf('[FAIL] fegeometry — %s\n', ME.message);
        fprintf('       identifier: %s\n\n', ME.identifier);
    end
end
