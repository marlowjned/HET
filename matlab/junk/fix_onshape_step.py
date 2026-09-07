#!/usr/bin/env python3
"""
fix_onshape_step.py

Restructures Onshape-exported STEP files so MATLAB's fegeometry can import them.

ROOT CAUSE: Onshape wraps every export in an assembly structure. MATLAB's STEP
reader requires a flat product structure (confirmed against MATLAB's own example
STEP files). The assembly wrapper entities cause the STEPGeometryMustBeSolid error.

WHAT THIS SCRIPT REMOVES:
  - NEXT_ASSEMBLY_USAGE_OCCURRENCE   (assembly instance links)
  - CONTEXT_DEPENDENT_SHAPE_REPRESENTATION  (instance placement contexts)
  - Compound REPRESENTATION_RELATIONSHIP + SHAPE_REPRESENTATION_RELATIONSHIP
    entities  (instance placement transformations)
  - ITEM_DEFINED_TRANSFORMATION      (placement matrices)
  - The root assembly PRODUCT_DEFINITION chain  (PD, PDS, SDR)

WHAT IT KEEPS:
  - All part-level PRODUCT → PDF → PD → PDS → SDR → SHAPE_REPRESENTATION chains
  - All SHAPE_REPRESENTATION_RELATIONSHIP entities (the SRR linking plain rep → ABSR)
  - All ADVANCED_BREP_SHAPE_REPRESENTATION and solid geometry entities

Usage:
  python fix_onshape_step.py input.step [output.step]
  python fix_onshape_step.py *.step          (batch: writes *_clean.step)
"""

import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Entity block parser
# ---------------------------------------------------------------------------

def parse_blocks(text):
    """Split STEP DATA section into a list of blocks.

    Each block is a dict:
        {'type': 'entity'|'other', 'lines': [str, ...], 'id': int|None}
    """
    lines = text.splitlines(keepends=True)
    blocks = []
    i = 0
    while i < len(lines):
        raw = lines[i]
        stripped = raw.strip()

        if not stripped.startswith('#'):
            blocks.append({'type': 'other', 'lines': [raw], 'id': None})
            i += 1
            continue

        # Start of an entity definition
        entity_lines = [raw]
        if stripped.endswith(';'):
            # Single-line entity
            eid = _leading_id(stripped)
            blocks.append({'type': 'entity', 'lines': entity_lines, 'id': eid})
            i += 1
        else:
            # Multi-line entity: accumulate until a line ending with ');' or ');'
            i += 1
            while i < len(lines):
                entity_lines.append(lines[i])
                s = lines[i].strip()
                if s.endswith(');') or s == ');':
                    break
                i += 1
            eid = _leading_id(entity_lines[0].strip())
            blocks.append({'type': 'entity', 'lines': entity_lines, 'id': eid})
            i += 1

    return blocks


def _leading_id(line):
    m = re.match(r'^#(\d+)\s*=', line)
    return int(m.group(1)) if m else None


def entity_type(block):
    """Return the STEP entity type string, or 'COMPOUND' for (#N=(...))."""
    first = block['lines'][0].strip()
    m = re.match(r'^#\d+\s*=\s*([A-Z_]+)\s*\(', first)
    if m:
        return m.group(1)
    if re.match(r'^#\d+\s*=\s*\(', first):
        return 'COMPOUND'
    return None


def all_refs(block):
    """Return all #N references that appear in the block."""
    text = ''.join(block['lines'])
    return [int(x) for x in re.findall(r'#(\d+)', text)]


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def fix_step(input_path, output_path=None):
    if output_path is None:
        p = Path(input_path)
        output_path = p.with_name(p.stem + '_clean.step')

    text = Path(input_path).read_text(errors='replace')
    blocks = parse_blocks(text)

    # Index entities by id
    by_id = {}
    for b in blocks:
        if b['type'] == 'entity' and b['id'] is not None:
            by_id[b['id']] = b

    print(f'\n=== {Path(input_path).name} ===')

    # ------------------------------------------------------------------
    # 1. Identify entity IDs to remove by type
    # ------------------------------------------------------------------
    remove = set()

    nauo_ids      = set()
    cdsr_ids      = set()
    idt_ids       = set()
    comp_srr_ids  = set()   # compound SRR/RR placement entities

    for eid, b in by_id.items():
        t = entity_type(b)
        if t == 'NEXT_ASSEMBLY_USAGE_OCCURRENCE':
            nauo_ids.add(eid)
        elif t == 'CONTEXT_DEPENDENT_SHAPE_REPRESENTATION':
            cdsr_ids.add(eid)
        elif t == 'ITEM_DEFINED_TRANSFORMATION':
            idt_ids.add(eid)
        elif t == 'COMPOUND':
            full = ''.join(b['lines'])
            if 'SHAPE_REPRESENTATION_RELATIONSHIP()' in full:
                comp_srr_ids.add(eid)

    remove |= nauo_ids | cdsr_ids | idt_ids | comp_srr_ids
    print(f'  Assembly instances  NAUO:{sorted(nauo_ids)}  CDSR:{sorted(cdsr_ids)}'
          f'  IDT:{sorted(idt_ids)}  compSRR:{sorted(comp_srr_ids)}')

    # ------------------------------------------------------------------
    # 2. Find root assembly PRODUCT_DEFINITION chain
    #    Root = the PD that is a NAUO *parent* but never a NAUO *child*.
    # ------------------------------------------------------------------

    # NAUO format: NEXT_ASSEMBLY_USAGE_OCCURRENCE('id','name','desc',
    #              #relating_pd, #related_pd, 'ref_des')
    # The second-to-last #ref is the relating (parent) PD.
    # The last #ref is the related (child) PD.
    nauo_parent_pds = set()
    nauo_child_pds  = set()
    for eid in nauo_ids:
        refs = all_refs(by_id[eid])
        # refs[0] = eid itself; last two non-eid refs are parent, child
        non_self = [r for r in refs if r != eid]
        if len(non_self) >= 2:
            nauo_parent_pds.add(non_self[-2])
            nauo_child_pds.add(non_self[-1])

    root_pds = nauo_parent_pds - nauo_child_pds
    print(f'  Root assembly PD(s): {sorted(root_pds)}')

    # For each root PD, trace: PD -> PDS -> SDR
    # PDS: PRODUCT_DEFINITION_SHAPE(name, desc, #pd)
    # SDR: SHAPE_DEFINITION_REPRESENTATION(#pds, #shape_rep)

    pds_for_pd = {}   # pd_id -> [pds_id, ...]
    sdr_for_pds = {}  # pds_id -> [sdr_id, ...]

    for eid, b in by_id.items():
        t = entity_type(b)
        if t == 'PRODUCT_DEFINITION_SHAPE':
            refs = all_refs(b)
            if refs:
                pd_ref = refs[-1]
                pds_for_pd.setdefault(pd_ref, []).append(eid)
        elif t == 'SHAPE_DEFINITION_REPRESENTATION':
            refs = all_refs(b)
            if len(refs) >= 2:
                pds_ref = refs[1]
                sdr_for_pds.setdefault(pds_ref, []).append(eid)

    root_pds_ids  = set()
    root_sdr_ids  = set()
    root_nauo_pds_ids = set()  # PDS for NAUO occurrences (#516 type)

    for pd in root_pds:
        for pds in pds_for_pd.get(pd, []):
            root_pds_ids.add(pds)
            for sdr in sdr_for_pds.get(pds, []):
                root_sdr_ids.add(sdr)

    # Also remove PDS entities that reference NAUO entities directly
    for eid, b in by_id.items():
        if entity_type(b) == 'PRODUCT_DEFINITION_SHAPE':
            refs = all_refs(b)
            if refs and refs[-1] in nauo_ids:
                root_nauo_pds_ids.add(eid)

    remove |= root_pds_ids | root_sdr_ids | root_nauo_pds_ids
    print(f'  Root assembly PDS:{sorted(root_pds_ids)}  SDR:{sorted(root_sdr_ids)}'
          f'  NAUO-PDS:{sorted(root_nauo_pds_ids)}')
    print(f'  Total entities to remove: {len(remove)}')

    # ------------------------------------------------------------------
    # 3. Write output, skipping removed entities
    # ------------------------------------------------------------------
    out_parts = []
    skipped = 0
    for b in blocks:
        if b['type'] == 'entity' and b['id'] in remove:
            skipped += 1
            continue
        out_parts.append(''.join(b['lines']))

    Path(output_path).write_text(''.join(out_parts), encoding='utf-8')
    print(f'  Removed {skipped} entities  ->  {output_path}')
    return str(output_path)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        sys.exit(1)

    # Expand simple glob patterns on Windows (where the shell doesn't)
    import glob
    paths = []
    for a in args:
        expanded = glob.glob(a)
        paths.extend(expanded if expanded else [a])

    for p in paths:
        if not Path(p).exists():
            print(f'[SKIP] not found: {p}')
            continue
        fix_step(p)
