// magnetostatics_toy.pro
//
// 3-D magnetostatics via the magnetic vector potential (Whitney edge
// elements), single-coil validation case. Solved domain is the whole air
// sphere (simply connected) with a=0 Dirichlet on the true outer boundary
// -- no explicit tree/cohomology gauge needed, since the null space of
// curl-curl on a simply connected domain with full Dirichlet closure is
// just grad(0). Compare the resulting on-axis B at the coil center against
// the ideal-solenoid estimate B = mu0*N*I/L printed by build_toy.py.

Group {
  // Tags pinned explicitly in build_toy.py (GetDP Region[] takes integer
  // tags only, not physical-name strings): 1=Air, 2=Coil, 3=OuterBnd
  Air      = Region[1];
  Coil     = Region[2];
  OuterBnd = Region[3];

  DomainC  = Region[{Coil}];
  DomainCC = Region[{Air}];
  Domain   = Region[{DomainC, DomainCC}];
}

Function {
  mu0 = 4*Pi*1e-7;
  nu[Air]  = 1/mu0;
  nu[Coil] = 1/mu0;

  Jmag = 5.638644e+05;  // A/m^2, from build_toy.py (turns*current/A_cross)
  Js[Coil] = (Jmag / Sqrt[X[]^2 + Y[]^2]) * Vector[-Y[], X[], 0];
}

Jacobian {
  { Name JVol;
    Case {
      { Region All; Jacobian Vol; }
    }
  }
}

Integration {
  { Name Int1;
    Case {
      { Type Gauss;
        Case {
          { GeoElement Tetrahedron; NumberOfPoints 4; }
          { GeoElement Triangle; NumberOfPoints 3; }
        }
      }
    }
  }
}

Constraint {
  { Name MVP_0;
    Case {
      { Region OuterBnd; Value 0; }
    }
  }
  // Tree-cotree gauge: curl-curl's null space is grad(phi) for any phi in
  // H1_0(Domain) (not just harmonic phi -- a=0 on OuterBnd alone does NOT
  // remove this, since Domain's interior is unconstrained), so without this
  // the system is singular. Pin "a" to zero on a spanning tree of Domain's
  // edges, rooted on OuterBnd (where the Dirichlet constraint above already
  // holds) -- this kills the gradient null space without touching the
  // physical (solenoidal) part of "a". See GetDP manual tutorials 5/7/10.
  { Name MVP_Gauge;
    Case {
      { Region Domain; SubRegion OuterBnd; Value 0.; }
    }
  }
}

FunctionSpace {
  { Name Hcurl_a; Type Form1;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_Edge;
        Support Domain; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef ae; EntityType EdgesOf; NameOfConstraint MVP_0; }
      { NameOfCoef ae; EntityType EdgesOfTreeIn; EntitySubType StartingOn;
        NameOfConstraint MVP_Gauge; }
    }
  }
}

Formulation {
  { Name Magnetostatics_a; Type FemEquation;
    Quantity {
      { Name a; Type Local; NameOfSpace Hcurl_a; }
    }
    Equation {
      Galerkin { [ nu[] * Dof{d a}, {d a} ];
        In Domain; Jacobian JVol; Integration Int1; }
      Galerkin { [ -Js[], {a} ];
        In DomainC; Jacobian JVol; Integration Int1; }
    }
  }
}

Resolution {
  { Name Res_a;
    System {
      { Name Sys_a; NameOfFormulation Magnetostatics_a; }
    }
    Operation {
      Generate[Sys_a]; Solve[Sys_a]; SaveSolution[Sys_a];
    }
  }
}

PostProcessing {
  { Name PostPro_a; NameOfFormulation Magnetostatics_a;
    Quantity {
      { Name b; Value { Local { [ {d a} ]; In Domain; Jacobian JVol; } } }
    }
  }
}

PostOperation {
  { Name Map_b; NameOfPostProcessing PostPro_a;
    Operation {
      Print[ b, OnElementsOf Domain, File "b_toy.pos" ];
      Print[ b, OnPoint {0,0,0}, Format Table, File "b_center.txt" ];
    }
  }
}
