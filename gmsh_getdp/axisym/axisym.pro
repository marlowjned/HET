// axisym.pro
//
// 2-D axisymmetric linear magnetostatics, a-formulation (a = a_phi, Form1P
// perpendicular edges) with the VolAxiSqu Jacobian -- the axisymmetric branch
// of GetDP's own templates/Lib_Magnetostatics_a_phi.pro. Gmsh plane:
// x = r, y = thrust axis. Regions and winding areas come from regions.pro,
// written by geometry.py.
//
// Excitation in ampere-turns, set per run: -setnumber NI_INNER 1 -setnumber
// NI_OUTER 0 gives the inner basis field, and vice versa. NI_OUTER is per
// outer coil; the two outer sheets carry +-NI_OUTER (see geometry.py). Sign
// convention: both bases give B_r > 0 (radially outward) in the channel,
// which calibrate.py checks. NI_TRIM drives the optional trim coil.

Include "regions.pro";

DefineConstant[ NI_INNER = 1, NI_OUTER = 0, NI_TRIM = 0, Map = 0 ];

mu0 = 4e-7 * Pi;

Group {
  Windings = Region[{WInner, WOutIn, WOutOut, WTrim}];
  Domain = Region[{Air, Iron, Ring, Windings}];
}

Function {
  nu[Region[{Air, Windings}]] = 1 / mu0;
  nu[Iron] = 1 / (mu0 * MUR_IRON);
  nu[Ring] = 1 / (mu0 * MUR_RING);
  js[WInner]  = Vector[0, 0, -NI_INNER / A_WIN];
  js[WOutIn]  = Vector[0, 0, -NI_OUTER / A_WOI];
  js[WOutOut] = Vector[0, 0,  NI_OUTER / A_WOO];
  // trim coil: same sense as the inner coil for NI_TRIM > 0
  js[WTrim]   = Vector[0, 0, -NI_TRIM / A_TRM];
}

Jacobian {
  { Name JVol; Case { { Region All; Jacobian VolAxiSqu; } } }
}

Integration {
  { Name Int; Case { { Type Gauss; Case { { GeoElement Triangle; NumberOfPoints 4; } } } } }
}

Constraint {
  { Name A0; Case { { Region Bnd; Value 0; } } }
}

FunctionSpace {
  { Name Hcurl_a; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge; Support Domain; Entity NodesOf[All]; }
    }
    Constraint { { NameOfCoef ae; EntityType NodesOf; NameOfConstraint A0; } }
  }
}

Formulation {
  { Name MS; Type FemEquation;
    Quantity { { Name a; Type Local; NameOfSpace Hcurl_a; } }
    Equation {
      Integral { [ nu[] * Dof{d a}, {d a} ]; In Domain; Jacobian JVol; Integration Int; }
      Integral { [ -js[], {a} ]; In Windings; Jacobian JVol; Integration Int; }
    }
  }
}

Resolution {
  { Name Res; System { { Name S; NameOfFormulation MS; } }
    Operation { Generate[S]; Solve[S]; SaveSolution[S]; } }
}

PostProcessing {
  { Name PP; NameOfFormulation MS;
    Quantity {
      { Name b;   Value { Local { [ {d a} ]; In Domain; Jacobian JVol; } } }
      // r * a_phi: its contours are the field lines
      { Name raz; Value { Local { [ CompZ[{a}] * X[] ]; In Domain; Jacobian JVol; } } }
    }
  }
}

// Channel sampling grid: r across the channel, y from 30 mm downstream of the
// exit (-67.846) to 3 mm behind the anode face (-1.826), 0.25 mm steps.
PostOperation {
  { Name Post; NameOfPostProcessing PP;
    Operation {
      Print[ b, OnGrid {$A, $B, 0} { {0.022120:0.031600:0.000948}, {-0.067846:-0.001826:0.00025}, {0} },
             File "b_channel.txt", Format SimpleTable ];
      Print[ b, OnElementsOf Region[{Iron, Ring}], Depth 0, File "b_iron.txt", Format SimpleTable ];
      If(Map)
        // element-centroid values (Depth 0): far cheaper than OnGrid, which
        // locates every grid point in the mesh one by one
        Print[ b, OnElementsOf Domain, Depth 0, File "b_map.txt", Format SimpleTable ];
        Print[ raz, OnElementsOf Domain, Depth 0, File "raz_map.txt", Format SimpleTable ];
      EndIf
    }
  }
}
