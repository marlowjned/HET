%% het_cad_diag.m  — temporary diagnostic, safe to delete after
clear; close all; clc;

base = fileparts(mfilename('fullpath'));

%% Test 1: fegeometry with analytic geometry (no STEP)
try
    gm = multicylinder(0.05, 0.1);
    gm = fegeometry(gm);
    fprintf('[PASS] Test 1: fegeometry(analytic) OK  (NumCells=%d)\n', gm.NumCells);
catch ME
    fprintf('[FAIL] Test 1: fegeometry(analytic) — %s\n', ME.message);
end

%% Test 2: importGeometry legacy API (single part, AP214)
try
    model = createpde(1);
    importGeometry(model, fullfile(base, 'BPL-700-part-chamber.step'));
    fprintf('[PASS] Test 2: importGeometry(chamber AP214) OK  (NumCells=%d)\n', model.Geometry.NumCells);
catch ME
    fprintf('[FAIL] Test 2: importGeometry(chamber AP214) — %s\n', ME.message);
end

%% Test 3: fegeometry on single part (AP214)
try
    gm = fegeometry(fullfile(base, 'BPL-700-part-chamber.step'));
    fprintf('[PASS] Test 3: fegeometry(chamber AP214) OK  (NumCells=%d)\n', gm.NumCells);
catch ME
    fprintf('[FAIL] Test 3: fegeometry(chamber AP214) — %s\n', ME.message);
    fprintf('       identifier: %s\n', ME.identifier);
    if ~isempty(ME.cause)
        fprintf('       cause: %s\n', ME.cause{1}.message);
    end
end

%% Test 4: fegeometry on single part (AP203)
try
    gm = fegeometry(fullfile(base, 'BPL-700-part-injector.step'));
    fprintf('[PASS] Test 4: fegeometry(injector AP214) OK  (NumCells=%d)\n', gm.NumCells);
catch ME
    fprintf('[FAIL] Test 4: fegeometry(injector AP214) — %s\n', ME.message);
    fprintf('       identifier: %s\n', ME.identifier);
end

%% Test 5: fegeometry on AP214 assembly (6 bodies)
try
    gm = fegeometry(fullfile(base, 'BPL-700-sims-assem.step'));
    fprintf('[PASS] Test 5: fegeometry(AP214 assembly) OK  (NumCells=%d)\n', gm.NumCells);
catch ME
    fprintf('[FAIL] Test 5: fegeometry(AP214 assembly) — %s\n', ME.message);
    fprintf('       identifier: %s\n', ME.identifier);
end

%% Test 6: fegeometry on AP203 assembly (13 bodies)
try
    gm = fegeometry(fullfile(base, 'BPL-700-sims-assem-203.step'));
    fprintf('[PASS] Test 6: fegeometry(AP203 assembly) OK  (NumCells=%d)\n', gm.NumCells);
catch ME
    fprintf('[FAIL] Test 6: fegeometry(AP203 assembly) — %s\n', ME.message);
    fprintf('       identifier: %s\n', ME.identifier);
end

%% Summary
fprintf('\nMATLAB version: %s\n', version);
fprintf('PDE Toolbox:    '); ver('pde');
