#!/usr/bin/env python3
"""Extract per-NAUO placement transforms from an Onshape assembly STEP file
to verify whether patterned instances (e.g. the 4 outer_solenoid copies)
actually carry distinct rotation/translation data, and print each with its
child part name so specific occurrences can be identified."""
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text(errors='replace')
lines = text.splitlines()

def parse_entities(lines):
    ents = {}
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        m = re.match(r'^#(\d+)\s*=\s*(\w+)\s*\((.*)$', s)
        if not m:
            m2 = re.match(r'^#(\d+)\s*=\s*(\(.*)$', s)   # compound entity: #N=(...)
            if m2:
                eid = int(m2.group(1))
                buf = m2.group(2)
                j = i
                while not buf.rstrip().endswith(');') and j < len(lines) - 1:
                    j += 1
                    buf += lines[j].strip()
                ents[eid] = ('COMPOUND', buf)
                i = j + 1
                continue
            i += 1
            continue
        eid = int(m.group(1))
        etype = m.group(2)
        buf = m.group(3)
        j = i
        while not buf.rstrip().endswith(');') and j < len(lines) - 1:
            j += 1
            buf += lines[j].strip()
        ents[eid] = (etype, buf)
        i = j + 1
    return ents

ents = parse_entities(lines)

def refs(s):
    return [int(x) for x in re.findall(r'#(\d+)', s)]

def floats(s):
    pat = r'[-+]?(?:\d+\.\d*|\.\d+)(?:E[-+]?\d+)?'
    return [float(x) for x in re.findall(pat, s)]

nauos  = {eid: buf for eid, (t, buf) in ents.items() if t == 'NEXT_ASSEMBLY_USAGE_OCCURRENCE'}
pdss   = {eid: buf for eid, (t, buf) in ents.items() if t == 'PRODUCT_DEFINITION_SHAPE'}
cdsrs  = {eid: buf for eid, (t, buf) in ents.items() if t == 'CONTEXT_DEPENDENT_SHAPE_REPRESENTATION'}
idts   = {eid: buf for eid, (t, buf) in ents.items() if t == 'ITEM_DEFINED_TRANSFORMATION'}
comps  = {eid: buf for eid, (t, buf) in ents.items() if t == 'COMPOUND'}
axis_placements = {eid: buf for eid, (t, buf) in ents.items() if t == 'AXIS2_PLACEMENT_3D'}
cart_pts = {eid: buf for eid, (t, buf) in ents.items() if t == 'CARTESIAN_POINT'}
dirs = {eid: buf for eid, (t, buf) in ents.items() if t == 'DIRECTION'}

print(f'Found {len(nauos)} NAUO, {len(pdss)} PDS, {len(cdsrs)} CDSR, {len(comps)} COMPOUND, {len(idts)} IDT')

def placement_full(axId):
    if axId not in axis_placements:
        return None
    r = refs(axis_placements[axId])
    loc  = r[0] if len(r) > 0 else None
    axis = r[1] if len(r) > 1 else None
    refd = r[2] if len(r) > 2 else None
    return {
        'loc':    floats(cart_pts.get(loc, ''))  if loc  else None,
        'axis':   floats(dirs.get(axis, ''))     if axis else None,
        'refdir': floats(dirs.get(refd, ''))     if refd else None,
    }

nauo_child_name = {}
for eid, buf in nauos.items():
    m = re.search(r"'([^']*)'\s*,\s*'([^']*)'", buf)
    nauo_child_name[eid] = m.group(2) if m else '?'

# PDS -> NAUO: PDS's last ref is the NAUO id (this PDS represents "the shape of this occurrence")
pds_to_nauo = {}
for eid, buf in pdss.items():
    r = refs(buf)
    if r and r[-1] in nauos:
        pds_to_nauo[eid] = r[-1]

# CDSR -> (compound_id, pds_id)
naou_to_idt = {}
for cdsr_eid, buf in cdsrs.items():
    r = refs(buf)
    compound_id = next((x for x in r if x in comps), None)
    pds_id      = next((x for x in r if x in pds_to_nauo), None)
    if compound_id is None or pds_id is None:
        continue
    nauo_id = pds_to_nauo[pds_id]
    idt_id  = next((x for x in refs(comps[compound_id]) if x in idts), None)
    if idt_id is not None:
        naou_to_idt[nauo_id] = idt_id

print(f'\nResolved {len(naou_to_idt)} NAUO -> IDT links\n')
for nauo_id, idt_id in naou_to_idt.items():
    r = refs(idts[idt_id])
    p2 = r[-1]
    pf = placement_full(p2)
    name = nauo_child_name.get(nauo_id, '?')
    print(f'  NAUO#{nauo_id} ({name}) -> IDT#{idt_id}: loc={pf["loc"]}  axis(Z)={pf["axis"]}  refdir(X)={pf["refdir"]}')
