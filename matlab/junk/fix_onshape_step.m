function fix_onshape_step(inputFile, outputFile)
% fix_onshape_step  Patch an Onshape-exported STEP file for MATLAB fegeometry.
%
% Problem: Onshape wraps every solid in an assembly structure. MATLAB's STEP
% reader finds SHAPE_DEFINITION_REPRESENTATION pointing to a plain SHAPE_
% REPRESENTATION (axis placements only), never reaching the solid geometry.
%
% Fix: for each SHAPE_DEFINITION_REPRESENTATION that points to a plain rep
% which has a SHAPE_REPRESENTATION_RELATIONSHIP leading to an ADVANCED_BREP,
% redirect it to point directly to the ADVANCED_BREP instead.
%
% Usage:
%   fix_onshape_step('in.step')            % writes in_fixed.step
%   fix_onshape_step('in.step','out.step')

    if nargin < 2
        [d, n, e] = fileparts(inputFile);
        outputFile = fullfile(d, [n '_fixed' e]);
    end

    text = fileread(inputFile);
    lines = splitlines(text);
    fprintf('Processing %s (%d lines)\n', inputFile, numel(lines));

    %% 1. Find ADVANCED_BREP_SHAPE_REPRESENTATION entity IDs
    absr_ids = [];
    for i = 1:numel(lines)
        ln = strtrim(lines{i});
        if contains(ln, '=ADVANCED_BREP_SHAPE_REPRESENTATION(')
            id = extractLeadingId(ln);
            if ~isnan(id), absr_ids(end+1) = id; end %#ok<AGROW>
        end
    end
    if isempty(absr_ids)
        error('No ADVANCED_BREP_SHAPE_REPRESENTATION found — not a solid STEP file?');
    end
    fprintf('  Solid bodies (ADVANCED_BREP): %s\n', num2str(absr_ids));

    %% 2. Build map: plain_shape_rep_id -> absr_id
    %     via standalone SHAPE_REPRESENTATION_RELATIONSHIP lines.
    %     Standalone SRRs look like:
    %       #N=SHAPE_REPRESENTATION_RELATIONSHIP('...','...',#rep1,#rep2);
    %     Compound SRRs (to be ignored) look like:
    %       SHAPE_REPRESENTATION_RELATIONSHIP()     <- inside a (#N=(...); block
    plainToAbsr = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    for i = 1:numel(lines)
        ln = strtrim(lines{i});
        % Must start with '#' and contain the full entity type name with args
        if ~startsWith(ln, '#'), continue; end
        if ~contains(ln, '=SHAPE_REPRESENTATION_RELATIONSHIP('), continue; end
        % Skip the empty compound form: SHAPE_REPRESENTATION_RELATIONSHIP()
        if contains(ln, 'SHAPE_REPRESENTATION_RELATIONSHIP()'), continue; end

        % Extract all #N references from this line; last two are rep1, rep2
        refs = regexp(ln, '#(\d+)', 'tokens');
        if numel(refs) < 3, continue; end   % need entity id + 2 args
        rep1 = int32(str2double(refs{end-1}{1}));
        rep2 = int32(str2double(refs{end}{1}));
        if ismember(double(rep2), absr_ids)
            plainToAbsr(rep1) = rep2;
        elseif ismember(double(rep1), absr_ids)
            plainToAbsr(rep2) = rep1;
        end
    end

    if plainToAbsr.Count == 0
        warning('fix_onshape_step:noSRR', ...
            'No SRR -> ADVANCED_BREP links found. Structure may already be flat.');
    end
    k = keys(plainToAbsr); v = values(plainToAbsr);
    for i = 1:numel(k)
        fprintf('  SRR: plain rep #%d -> ABSR #%d\n', k{i}, v{i});
    end

    %% 3. Update SHAPE_DEFINITION_REPRESENTATION entries
    %     Format: #N=SHAPE_DEFINITION_REPRESENTATION(#prod_def,#shape_rep);
    nUpdated = 0;
    for i = 1:numel(lines)
        ln = strtrim(lines{i});
        if ~startsWith(ln, '#'), continue; end
        if ~contains(ln, '=SHAPE_DEFINITION_REPRESENTATION('), continue; end

        refs = regexp(ln, '#(\d+)', 'tokens');
        if numel(refs) < 3, continue; end  % need entity id + prod_def + shape_rep

        shapeRep = int32(str2double(refs{3}{1}));
        if ~isKey(plainToAbsr, shapeRep), continue; end

        newRep = plainToAbsr(shapeRep);
        oldRef = sprintf(',#%d)', double(shapeRep));
        newRef = sprintf(',#%d)', double(newRep));
        lines{i} = strrep(ln, oldRef, newRef);

        fprintf('  SDR #%s: rep #%d -> #%d (ADVANCED_BREP)\n', ...
            refs{1}{1}, shapeRep, newRep);
        nUpdated = nUpdated + 1;
    end

    %% 3b. Also redirect any remaining SDRs that STILL don't point to an ABSR.
    %      These are assembly-wrapper SDRs (e.g. the root 'BPL-700 Assembly'
    %      product) that have no SRR chain to follow.  MATLAB traverses from
    %      the root product and errors if its SDR doesn't lead to a solid, so
    %      we redirect them to the first available ABSR.
    firstAbsr = int32(absr_ids(1));
    for i = 1:numel(lines)
        ln = strtrim(lines{i});
        if ~startsWith(ln, '#'), continue; end
        if ~contains(ln, '=SHAPE_DEFINITION_REPRESENTATION('), continue; end

        refs = regexp(ln, '#(\d+)', 'tokens');
        if numel(refs) < 3, continue; end

        shapeRep = int32(str2double(refs{3}{1}));
        if ismember(double(shapeRep), absr_ids), continue; end  % already solid

        oldRef = sprintf(',#%d)', double(shapeRep));
        newRef = sprintf(',#%d)', double(firstAbsr));
        if strcmp(oldRef, newRef), continue; end
        lines{i} = strrep(ln, oldRef, newRef);
        fprintf('  SDR #%s: assembly wrapper #%d -> #%d (first ABSR)\n', ...
            refs{1}{1}, shapeRep, firstAbsr);
        nUpdated = nUpdated + 1;
    end

    if nUpdated == 0
        warning('fix_onshape_step:noUpdate', ...
            'No SHAPE_DEFINITION_REPRESENTATION entries updated. Check file structure.');
    else
        fprintf('  %d total SDR(s) updated.\n', nUpdated);
    end

    %% 4. Write output
    newText = strjoin(lines, newline);
    fid = fopen(outputFile, 'w');
    fwrite(fid, newText, 'char');
    fclose(fid);
    fprintf('  Saved: %s\n\n', outputFile);
end

function id = extractLeadingId(line)
    tok = regexp(line, '^#(\d+)\s*=', 'tokens', 'once');
    if isempty(tok), id = NaN; else, id = str2double(tok{1}); end
end
