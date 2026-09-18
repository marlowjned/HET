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
//
// Iron is a NONLINEAR material (generic soft-steel B-H curve, see
// generate_regions.py) -- this is why the Formulation below splits the
// single "nu[] * Dof{d a}" term into a linear part (Air+Windings) and a
// Picard-iteration nonlinear part (Iron, nu[{d a}]). Both the split and
// the per-stage residual loop are copied from GetDP's own bundled
// tools/getdp-3.5.0-Windows64/templates/Lib_Magnetostatics_a_phi.pro
// (Formulation "Magnetostatics_a" and its Resolution block) -- the
// standard idiom for a nonlinear a-formulation solve, not invented here.
//
// Plain Picard from a cold (zero-field) start DIVERGES on this curve --
// tried first, confirmed empirically: residual dropped 613->96 over 2
// iterations, then exploded to 1.7e5 on the 3rd. Root cause: this B-H
// curve's slope gets extremely steep near saturation (H jumps from
// 3.1e5 to 7.6e5 A/m over the curve's last few tenths of a Tesla), and
// the linear-material solve's own peak (~2.16 T, see gmsh_getdp/README.md)
// starts the very first nonlinear iteration right at that steep region --
// too large a jump for unrelaxed Picard to track. Fix: LOAD-STEPPING
// (continuation) -- ramp the excitation from a small fraction of full
// current up to 100% over several stages ($IFrac below, referenced in
// Js[] in assembly_regions_generated.pro), Picard-converging fully at
// each stage before stepping up. Each stage starts from the PREVIOUS
// stage's already-converged field (no InitSolution between stages), so
// every individual Picard step only has to track a small change instead
// of jumping from zero straight into saturation.

Include "assembly_regions_generated.pro";

DefineConstant[
  NL_tol_abs = 1e-6,   // absolute tolerance on residual for the nonlinear iron iteration
  NL_tol_rel = 1e-6,   // relative tolerance on residual for the nonlinear iron iteration
  NL_iter_max = 20     // maximum nonlinear (Picard) iterations PER load-step
];

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
      // Linear part: Air + Windings, nu[] is a plain per-region constant
      // (no argument needed).
      Galerkin { [ nu[] * Dof{d a}, {d a} ];
        In Region[{Air, Windings}]; Jacobian JVol; Integration Int1; }
      // Nonlinear part: Iron, nu[{d a}] evaluated at the CURRENT solution's
      // field (Picard fixed-point form -- simpler than full Newton, which
      // would additionally need a dhdb[]/dnudb2[] tangent term; see
      // generate_regions.py's docstring for why this is the deliberate v1
      // choice). Each resolution iteration below regenerates this term
      // using the previous iteration's {d a} as the linearization point.
      Galerkin { [ nu[{d a}] * Dof{d a}, {d a} ];
        In Iron; Jacobian JVol; Integration Int1; }
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
      // Load-stepping (continuation): ramp $IFrac 1/N_STEPS -> 1, Picard-
      // converging fully at each stage (Generate/Solve/GetResidual/While,
      // same pattern as GetDP's own Lib_Magnetostatics_a_phi.pro template)
      // before increasing the excitation further. See the header comment
      // above for why a single cold-start Picard solve diverges on this
      // B-H curve. InitSolution is called ONCE, before the first (smallest)
      // stage -- every later stage continues from the previous stage's
      // converged field, not from zero.
      InitSolution[Sys_a];
      For n In {1:8}
        Evaluate[ $IFrac = n / 8.0 ];
        Generate[Sys_a]; Solve[Sys_a];
        Generate[Sys_a]; GetResidual[Sys_a, $res0];
        Evaluate[ $res = $res0, $iter = 0 ];
        Print[{n, $iter, $res, $res / $res0},
          Format "IFrac step %g, NL iter %03g: residual abs %14.12e rel %14.12e"];
        While[$res > NL_tol_abs && $res / $res0 > NL_tol_rel &&
              $res / $res0 <= 1 && $iter < NL_iter_max]{
          Solve[Sys_a]; Generate[Sys_a]; GetResidual[Sys_a, $res];
          Evaluate[ $iter = $iter + 1 ];
          Print[{n, $iter, $res, $res / $res0},
            Format "IFrac step %g, NL iter %03g: residual abs %14.12e rel %14.12e"];
        }
      EndFor
      SaveSolution[Sys_a];
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
