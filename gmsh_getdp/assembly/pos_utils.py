"""
pos_utils.py

Shared parser for GetDP's Gmsh .pos vector-view output (used by
post_process.py, cross_section_viewer.py, and export_vtk.py).
"""
import numpy as np


def parse_pos_elements(path):
    """Returns (nodes[N,4,3], B[N,3]) -- one row per element: its 4 corner
    node coordinates, and its (P0 -- constant per element, see
    HOW_IT_WORKS.md sec. 9) field value."""
    nodes_list = []
    Bs = []
    with open(path) as f:
        for line in f:
            if not line.startswith("VS("):
                continue
            coord_str, val_str = line.split("){", 1)
            coords = [float(x) for x in coord_str[3:].split(",")]
            nodes = np.array(coords).reshape(4, 3)
            vals = [float(x) for x in val_str.rstrip(";\n").rstrip("}").split(",")]
            b = vals[0:3]  # same value repeated at all 4 nodes -- P0 field
            nodes_list.append(nodes)
            Bs.append(b)
    return np.array(nodes_list), np.array(Bs)


def parse_pos(path):
    """Returns (centroids[N,3], B[N,3]) -- one row per element."""
    nodes, B = parse_pos_elements(path)
    centroids = nodes.mean(axis=1)
    return centroids, B
