// magnetostatics_assembly.pro
//
// 3-D magnetostatics of the real BPL-700 assembly (iron + 5 windings + air),
// via the magnetic vector potential (Whitney edge elements). Same validated
// formulation as ../toy_test/magnetostatics_toy.pro (which matched an exact
// analytic on-axis field to 1.4%) -- see that file's comments for why the
// tree-cotree gauge below is required even though the domain is simply
// connected. Region tags, materials, and winding current sources are in the
// auto-generated assembly_regions_generated.pro (run build_assembly.py to
// regenerate after any geometry change).

Include "assembly_regions_generated.pro";

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
  // Tree-cotree gauge -- see magnetostatics_toy.pro for the full
  // explanation of why this is required.
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
      { Name b_iron; Value { Local { [ {d a} ]; In Iron; Jacobian JVol; } } }
    }
  }
}

PostOperation {
  { Name Map_b; NameOfPostProcessing PostPro_a;
    Operation {
      Print[ b, OnElementsOf Domain, File "b_assembly.pos" ];
      Print[ b_iron, OnElementsOf Iron, File "b_iron.pos" ];
    }
  }
}
