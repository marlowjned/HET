// estatic.pro: electrostatics of dielectrics + conductors, conductors either at an imposed potential or
// floating (imposed net charge, potential solved for). 2-D axisymmetric (x = r, y = z, axis x = 0) or 3-D.
// Scalar potential v with one global unknown per conductor (BF_GroupOfNodes), as GetDP's
// templates/Lib_Electrostatics_v.pro. Mesh in metres. Physical tags are set by mesh_rz.py / mesh_3d.py.
//
//   getdp estatic.pro -msh m.msh -solve ES -pos Post -setnumber V_anode 300 -setnumber Iron_floating 1
// Output files land next to the .pro file: run a copy of it in each work directory.

DefineConstant[
  Flag_Axi = 1,          // 1: 2-D axisymmetric; 0: 3-D
  eps_r_bn = 4.0,        // PLACEHOLDER relative permittivity of the BN
  V_anode = 1.0,         // [V]
  V_chamber = 0.0,       // [V] chamber wall (far boundary)
  Iron_floating = 1,     // 1: iron (with coils) floating at net charge Q_iron; 0: iron held at V_iron
  V_iron = 0.0,          // [V]
  Q_iron = 0.0,          // [C]
  Has_cathode = 0,       // 3-D only: a cathode conductor (tag 14)
  V_cathode = 0.0,
  P0x = 0, P0y = 0, P0z = 0, P1x = 1, P1y = 0, P1z = 0, P2x = 0, P2y = 0, P2z = 1, N1 = 10, N2 = 10  // Plane output
];

// GetDP's axisymmetric Jacobian integrates per radian: global charges and energies are scaled by 2 pi below,
// and an imposed charge is divided by it, so Q and energy are always totals (checked by verify_spheres.py)
If (Flag_Axi)
  Q_fac = 2 * Pi;
Else
  Q_fac = 1;
EndIf

Group {
  Vacuum = Region[1];
  BN = Region[2];
  Anode = Region[11];    // conductor surfaces
  Iron = Region[12];
  Chamber = Region[13];
  Cathode = Region[14];

  Vol_Ele = Region[{Vacuum, BN}];
  If (Has_cathode)
    Sur_C_Ele = Region[{Anode, Iron, Chamber, Cathode}];
  Else
    Sur_C_Ele = Region[{Anode, Iron, Chamber}];
  EndIf
}

Function {
  eps0 = 8.8541878128e-12;
  epsr[Vacuum] = 1.0;
  epsr[BN] = eps_r_bn;
}

Constraint {
  { Name GlobalPotential;
    Case {
      { Region Anode; Value V_anode; }
      { Region Chamber; Value V_chamber; }
      If (!Iron_floating)
        { Region Iron; Value V_iron; }
      EndIf
      If (Has_cathode)
        { Region Cathode; Value V_cathode; }
      EndIf
    }
  }
  { Name GlobalCharge;
    Case {
      If (Iron_floating)
        { Region Iron; Value Q_iron / Q_fac; }
      EndIf
    }
  }
  { Name NodalPotential; Case { } }
}

Jacobian {
  { Name Vol;
    Case {
      If (Flag_Axi)
        { Region All; Jacobian VolAxiSqu; }
      Else
        { Region All; Jacobian Vol; }
      EndIf
    }
  }
}

Integration {
  { Name Int;
    Case { { Type Gauss;
      Case {
        { GeoElement Point; NumberOfPoints 1; }
        { GeoElement Line; NumberOfPoints 3; }
        { GeoElement Triangle; NumberOfPoints 3; }
        { GeoElement Tetrahedron; NumberOfPoints 4; }
      }
    } }
  }
}

FunctionSpace {
  { Name Hgrad_v; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Vol_Ele; Entity NodesOf[ All, Not Sur_C_Ele ]; }
      { Name sc; NameOfCoef vc; Function BF_GroupOfNodes;
        Support Vol_Ele; Entity GroupsOfNodesOf[ Sur_C_Ele ]; }
    }
    GlobalQuantity {
      { Name GlobalPotential; Type AliasOf; NameOfCoef vc; }
      { Name GlobalCharge; Type AssociatedWith; NameOfCoef vc; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint NodalPotential; }
      { NameOfCoef GlobalPotential; EntityType GroupsOfNodesOf; NameOfConstraint GlobalPotential; }
      { NameOfCoef GlobalCharge; EntityType GroupsOfNodesOf; NameOfConstraint GlobalCharge; }
    }
  }
}

Formulation {
  { Name ES; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v; }
      { Name Q; Type Global; NameOfSpace Hgrad_v [GlobalCharge]; }
      { Name V; Type Global; NameOfSpace Hgrad_v [GlobalPotential]; }
    }
    Equation {
      Integral { [ epsr[] * eps0 * Dof{d v}, {d v} ]; In Vol_Ele; Jacobian Vol; Integration Int; }
      GlobalTerm { [ -Dof{Q}, {V} ]; In Sur_C_Ele; }
    }
  }
}

Resolution {
  { Name ES;
    System { { Name A; NameOfFormulation ES; } }
    Operation { Generate[A]; Solve[A]; SaveSolution[A]; }
  }
}

PostProcessing {
  { Name ES; NameOfFormulation ES;
    Quantity {
      { Name v; Value { Term { [ {v} ]; In Vol_Ele; Jacobian Vol; } } }
      { Name e; Value { Term { [ -{d v} ]; In Vol_Ele; Jacobian Vol; } } }
      { Name Q; Value { Term { [ Q_fac * {Q} ]; In Sur_C_Ele; } } }
      { Name V; Value { Term { [ {V} ]; In Sur_C_Ele; } } }
      { Name energy; Value { Integral { [ Q_fac * 0.5 * eps0 * epsr[] * SquNorm[{d v}] ];
                                        In Vol_Ele; Jacobian Vol; Integration Int; } } }
    }
  }
}

PostOperation {
  { Name Post; NameOfPostProcessing ES;
    Operation {
      Print[ v, OnElementsOf Vol_Ele, Depth 1, File "v_elem.txt", Format SimpleTable ];
      Print[ Q, OnRegion Anode, File "Q_anode.txt", Format Table ];
      Print[ V, OnRegion Anode, File "V_anode.txt", Format Table ];
      Print[ Q, OnRegion Iron, File "Q_iron.txt", Format Table ];
      Print[ V, OnRegion Iron, File "V_iron.txt", Format Table ];
      Print[ Q, OnRegion Chamber, File "Q_chamber.txt", Format Table ];
      If (Has_cathode)
        Print[ Q, OnRegion Cathode, File "Q_cathode.txt", Format Table ];
        Print[ V, OnRegion Cathode, File "V_cathode.txt", Format Table ];
      EndIf
      Print[ energy[Vol_Ele], OnGlobal, File "energy.txt", Format Table ];
    }
  }
  { Name Globals; NameOfPostProcessing ES;  // conductor charges and potentials only (3-D: no volume export)
    Operation {
      Print[ Q, OnRegion Anode, File "Q_anode.txt", Format Table ];
      Print[ V, OnRegion Anode, File "V_anode.txt", Format Table ];
      Print[ Q, OnRegion Iron, File "Q_iron.txt", Format Table ];
      Print[ V, OnRegion Iron, File "V_iron.txt", Format Table ];
      Print[ Q, OnRegion Chamber, File "Q_chamber.txt", Format Table ];
      If (Has_cathode)
        Print[ Q, OnRegion Cathode, File "Q_cathode.txt", Format Table ];
        Print[ V, OnRegion Cathode, File "V_cathode.txt", Format Table ];
      EndIf
    }
  }
  { Name Plane; NameOfPostProcessing ES;    // potential on a regular grid in a plane (origin P0, corners P1, P2)
    Operation {
      Print[ v, OnPlane { {P0x, P0y, P0z} {P1x, P1y, P1z} {P2x, P2y, P2z} } {N1, N2},
             File "plane.txt", Format SimpleTable ];
    }
  }
  { Name Maps; NameOfPostProcessing ES;     // for viewing in Gmsh
    Operation {
      Print[ v, OnElementsOf Vol_Ele, File "v.pos" ];
      Print[ e, OnElementsOf Vol_Ele, File "e.pos" ];
    }
  }
}
