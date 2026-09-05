# het_solenoid_bfield.m

Preliminary 3-D magnetostatic simulation of a Hall Effect Thruster electromagnet
layout: one center-pole solenoid plus four outer-pole solenoids spaced around
the discharge channel opening, driven at opposite polarity.

Builds the coil geometry, solves for the magnetic field with MATLAB's PDE
Toolbox, and plots the resulting field lines across the channel opening, plus
an interactive slider viewer for sweeping the field cross-section by azimuthal
angle. Used to rough-size the solenoids and confirm the field-generation
pipeline works before adding the thruster body/core material in a later pass.

Requires MATLAB R2023b+ with the Partial Differential Equation Toolbox.

## Future features / things to consider

- Full thruster material included
- Force acting on solenoids
- E field sim
- Trace electron trajectory, quantify stability
- Coil thermals?

