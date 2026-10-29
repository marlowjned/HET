// mag_rz.pro: 2-D axisymmetric linear magnetostatics, a-formulation (a = a_phi, Form1P perpendicular edges,
// VolAxiSqu Jacobian), as GetDP's templates/Lib_Magnetostatics_a_phi.pro and the HET 2-D magnet models.
// Gmsh plane: x = r, y = z. Mesh in metres; tags from run_mag.py.
//   getdp mag_rz.pro -msh m.msh -solve MS -pos Post -setnumber I_inner 1 -setnumber I_outer 0.79375

DefineConstant[
  I_inner = 1.0,          // [A] inner coil
  I_outer = 0.79375,      // [A] outer coils (q * 230/184 per amp of inner current)
  N_inner = 230, N_outer = 184,
  A_inner = 1, A_sheet = 1,        // [m^2] r-z areas of the coil regions (set by run_mag.py)
  mu_iron = 1500, mu_posts = 225.9
];

Group {
  Air = Region[1];
  Iron = Region[2];
  Posts = Region[3];
  CoilIn = Region[4];
  SheetIn = Region[5];
  SheetOut = Region[6];
  Bnd = Region[10];
  Windings = Region[{CoilIn, SheetIn, SheetOut}];
  Domain = Region[{Air, Iron, Posts, Windings}];
}

Function {
  mu0 = 4e-7 * Pi;
  nu[Region[{Air, Windings}]] = 1 / mu0;
  nu[Iron] = 1 / (mu0 * mu_iron);
  nu[Posts] = 1 / (mu0 * mu_posts);
  // J_theta = sign * turns * I / area; the out-of-plane axis of the gmsh plane (x = r, y = z) is -theta.
  // Inner coil and inner sheet +theta, outer sheet -theta (the 2-D smear of the outer coils).
  js[CoilIn] = Vector[0, 0, -N_inner * I_inner / A_inner];
  js[SheetIn] = Vector[0, 0, -N_outer * I_outer / A_sheet];
  js[SheetOut] = Vector[0, 0, N_outer * I_outer / A_sheet];
}

Jacobian { { Name JVol; Case { { Region All; Jacobian VolAxiSqu; } } } }

Integration {
  { Name Int; Case { { Type Gauss; Case { { GeoElement Triangle; NumberOfPoints 4; }
                                           { GeoElement Line; NumberOfPoints 4; } } } } }
}

Constraint { { Name A0; Case { { Region Bnd; Value 0; } } } }

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
  { Name MS; System { { Name S; NameOfFormulation MS; } }
    Operation { Generate[S]; Solve[S]; SaveSolution[S]; } }
}

PostProcessing {
  { Name PP; NameOfFormulation MS;
    Quantity {
      { Name b; Value { Local { [ {d a} ]; In Domain; Jacobian JVol; } } }
    }
  }
}

PostOperation {
  { Name Post; NameOfPostProcessing PP;
    Operation { Print[ b, OnElementsOf Domain, Depth 1, File "b_elem.txt", Format SimpleTable ]; }
  }
}
