%% het_solenoid_bfield.m
%
% Preliminary 3-D magnetostatic check for a HET electromagnet layout:
% 1 center-pole solenoid + N outer-pole solenoids around the discharge
% opening. No thruster body/core material yet (relative permeability = 1
% everywhere), so |B| here is a lower bound / topology check, not final.
% Center and outer coils are driven at opposite polarity so field lines
% arc across the annular gap instead of just adding axially.
%
% Requires MATLAB R2023b+ with PDE Toolbox (femodel/CellLoad/FaceBC).

clear; close all; clc;

%% 1. PARAMETERS (inches, converted to meters below -- SI needed for the solve)
in2m = 0.0254;

geom.OD_opening_in    = 4.0;
geom.ID_opening_in    = 3.0;
geom.depth_opening_in = 2.0;

centerCoil.OD_in     = 2.5;
centerCoil.ID_in     = 1.0;
centerCoil.height_in = 2.0;
centerCoil.turns     = 300;
centerCoil.current_A = 5.0;
centerCoil.sign      = +1;

outerCoil.n             = 4;
outerCoil.OD_in         = 1.25;
outerCoil.ID_in         = 0.75;
outerCoil.height_in     = 2.0;
outerCoil.clearance_in  = 0.25;  % radial gap from opening OD to coil ID
outerCoil.turns         = 200;
outerCoil.current_A     = 5.0;
outerCoil.sign          = -1;    % opposite polarity vs. center coil

airDomainFactor = 8;   % air sphere radius = this x max assembly extent

meshHmaxFactor  = 1;   % global element size = airRadius * this
meshHfineFactor = 1;   % coil-region element size = min coil wall * this

nGrid = 25;            % interpolation grid points per axis
nSeedR = 3;             % radial seed points across the channel gap
nSeedTheta = 12;        % azimuthal seed points

%% 2. DERIVED GEOMETRY (meters)
OD_opening    = geom.OD_opening_in * in2m;
ID_opening    = geom.ID_opening_in * in2m;
depth_opening = geom.depth_opening_in * in2m;

centerCoil.OD     = centerCoil.OD_in * in2m;
centerCoil.ID     = centerCoil.ID_in * in2m;
centerCoil.height = centerCoil.height_in * in2m;

outerCoil.OD        = outerCoil.OD_in * in2m;
outerCoil.ID        = outerCoil.ID_in * in2m;
outerCoil.height    = outerCoil.height_in * in2m;
outerCoil.clearance = outerCoil.clearance_in * in2m;

assert(centerCoil.OD < ID_opening, ...
    'Center coil OD must fit inside the inner opening diameter.');
outerCoil.R_centers = OD_opening/2 + outerCoil.clearance + outerCoil.OD/2;

centerCoil.A_cross = pi * ((centerCoil.OD/2)^2 - (centerCoil.ID/2)^2);
outerCoil.A_cross  = pi * ((outerCoil.OD/2)^2  - (outerCoil.ID/2)^2);

% Bulk equivalent current density for N wound turns carrying current I
centerCoil.J = centerCoil.sign * (centerCoil.turns * centerCoil.current_A) / centerCoil.A_cross;
outerCoil.J  = outerCoil.sign  * (outerCoil.turns  * outerCoil.current_A)  / outerCoil.A_cross;

assemblyExtentR = outerCoil.R_centers + outerCoil.OD/2;
assemblyExtentZ = max(centerCoil.height, outerCoil.height) / 2;
airRadius = airDomainFactor * max(assemblyExtentR, assemblyExtentZ);

fprintf('=== HET Solenoid Preliminary Sizing ===\n');
fprintf('Opening OD / ID / depth   : %.2f / %.2f / %.2f in\n', ...
    geom.OD_opening_in, geom.ID_opening_in, geom.depth_opening_in);
fprintf('Center coil OD/ID/height  : %.2f / %.2f / %.2f in, NI = %d A-turns, J = %.3g A/m^2\n', ...
    centerCoil.OD_in, centerCoil.ID_in, centerCoil.height_in, ...
    centerCoil.turns*centerCoil.current_A, centerCoil.J);
fprintf('Outer coil OD/ID/height   : %.2f / %.2f / %.2f in, NI = %d A-turns, J = %.3g A/m^2 (x%d, opposite polarity)\n', ...
    outerCoil.OD_in, outerCoil.ID_in, outerCoil.height_in, ...
    outerCoil.turns*outerCoil.current_A, outerCoil.J, outerCoil.n);
fprintf('Outer coil placement radius: %.2f in\n', outerCoil.R_centers/in2m);
fprintf('Air domain radius          : %.2f in\n\n', airRadius/in2m);

% Ideal-solenoid B ~ mu0*N*I/L, as a sanity check against the FEM result
mu0 = 4*pi*1e-7;
B_ideal_center = mu0 * centerCoil.turns * centerCoil.current_A / centerCoil.height;
B_ideal_outer  = mu0 * outerCoil.turns  * outerCoil.current_A  / outerCoil.height;
fprintf('Ideal on-axis |B| estimate : center %.4g mT, outer %.4g mT\n\n', ...
    B_ideal_center*1e3, B_ideal_outer*1e3);

%% 3. GEOMETRY
gm = multisphere(airRadius);

centerGm = multicylinder([centerCoil.ID/2, centerCoil.OD/2], centerCoil.height, ...
    Void=[true false]);
centerGm = translate(centerGm, [0, 0, -centerCoil.height/2]);
gm = addCell(gm, centerGm);
centerCellIdx = gm.NumCells;

outerCellIdx = zeros(1, outerCoil.n);
outerCoil.x0 = zeros(1, outerCoil.n);
outerCoil.y0 = zeros(1, outerCoil.n);
for k = 1:outerCoil.n
    theta_k = 2*pi*(k-1)/outerCoil.n;
    x0 = outerCoil.R_centers * cos(theta_k);
    y0 = outerCoil.R_centers * sin(theta_k);
    outerCoil.x0(k) = x0;
    outerCoil.y0(k) = y0;

    coilGm = multicylinder([outerCoil.ID/2, outerCoil.OD/2], outerCoil.height, ...
        Void=[true false]);
    coilGm = translate(coilGm, [x0, y0, -outerCoil.height/2]);
    gm = addCell(gm, coilGm);
    outerCellIdx(k) = gm.NumCells;
end

figure('Name', 'HET solenoid layout (geometry check)');
pdegplot(gm, FaceAlpha=0.2, CellLabels="on");
title('Geometry check: 1 center + N outer solenoids inside air domain');
drawnow;

%% 4. MODEL, MATERIAL, EXCITATION, BOUNDARY CONDITIONS
model = femodel(AnalysisType="magnetostatic", Geometry=gm);
model.VacuumPermeability = mu0;
model.MaterialProperties = materialProperties(RelativePermeability=1);

model.CellLoad(centerCellIdx) = cellLoad( ...
    CurrentDensity=@(region, state) windingCurrent(region, state, 0, 0, centerCoil.J));

for k = 1:outerCoil.n
    x0 = outerCoil.x0(k);
    y0 = outerCoil.y0(k);
    model.CellLoad(outerCellIdx(k)) = cellLoad( ...
        CurrentDensity=@(region, state) windingCurrent(region, state, x0, y0, outerCoil.J));
end

% Zero vector potential on the true outer air boundary (faces of Cell 1
% not shared with any coil cell)
allAirFaces  = cellFaces(model.Geometry, 1);
coilCellRange = [centerCellIdx, outerCellIdx];
coilFaces    = cellFaces(model.Geometry, coilCellRange);
outerBoundaryFaces = setdiff(allAirFaces, coilFaces);
model.FaceBC(outerBoundaryFaces) = faceBC(MagneticPotential=[0; 0; 0]);

%% 5. MESH + SOLVE
meshHmax  = airRadius * meshHmaxFactor;
minCoilWall = min(centerCoil.OD - centerCoil.ID, outerCoil.OD - outerCoil.ID) / 2;
meshHfine = minCoilWall * meshHfineFactor;

model = generateMesh(model, Hmax=meshHmax, Hface={coilFaces, meshHfine});

fprintf('Solving magnetostatic model (%d elements)...\n', size(model.Mesh.Elements, 2));
R = solve(model);
fprintf('Solve complete.\n\n');

%% 6. POST-PROCESSING
Bmag = sqrt(R.MagneticFluxDensity.Bx.^2 + ...
            R.MagneticFluxDensity.By.^2 + ...
            R.MagneticFluxDensity.Bz.^2);

fprintf('Peak |B| in solved domain : %.4g mT\n', max(Bmag)*1e3);

figure('Name', 'HET solenoid |B| on coil regions');
coilElem = findElements(R.Mesh, "region", Cell=coilCellRange);
pdeplot3D(R.Mesh.Nodes, R.Mesh.Elements(:, coilElem), ColorMapData=Bmag);
title('|B| on coil volumes');

% Grid kept smaller than the full air domain to keep interpolation cheap
gridExtentR = outerCoil.R_centers + outerCoil.OD;
gridExtentZ = max(centerCoil.height, outerCoil.height);
xg = linspace(-gridExtentR, gridExtentR, nGrid);
yg = linspace(-gridExtentR, gridExtentR, nGrid);
zg = linspace(-gridExtentZ, gridExtentZ, nGrid);
[X, Y, Z] = meshgrid(xg, yg, zg);

BI = R.interpolateMagneticFlux(X, Y, Z);
Bx = reshape(BI.Bx, size(X));
By = reshape(BI.By, size(Y));
Bz = reshape(BI.Bz, size(Z));

figure('Name', 'HET solenoid field lines');
pdegplot(gm, FaceAlpha=0.05);
hold on;

seedR = linspace(ID_opening/2, OD_opening/2, nSeedR);
seedTheta = linspace(0, 2*pi, nSeedTheta + 1);
seedTheta(end) = [];
[SR, STh] = meshgrid(seedR, seedTheta);
sx = SR(:) .* cos(STh(:));
sy = SR(:) .* sin(STh(:));
sz = zeros(size(sx));

BmagGrid = sqrt(Bx.^2 + By.^2 + Bz.^2);

% Trace each seed both along +B and -B so the full line through it is drawn
vertsFwd = stream3(X, Y, Z, Bx, By, Bz, sx, sy, sz);
vertsRev = stream3(X, Y, Z, -Bx, -By, -Bz, sx, sy, sz);
plotColoredStreamlines(vertsFwd, X, Y, Z, BmagGrid);
plotColoredStreamlines(vertsRev, X, Y, Z, BmagGrid);
plot3(sx, sy, sz, 'k.', MarkerSize=10);

% Scale pinned to the ideal-solenoid estimate, not raw min/max, since |B|
% spikes non-physically at the smeared current-density source regions
cb = colorbar;
cb.Label.String = '|B| (T)';
clim([0, 2*max(B_ideal_center, B_ideal_outer)]);

xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('B field lines colored by local |B|, seeded across the channel opening');
axis equal; view(35, 20);
hold off;

fprintf('Done. Inspect the geometry-check figure first to confirm coil placement,\n');
fprintf('then the |B| and field-line figures.\n');

%% 7. INTERACTIVE CROSS-SECTION VIEWER
% Slider sweeps the azimuthal cut angle theta (0-180 deg); r spans
% [-gridExtentR, gridExtentR] so both half-planes (theta and theta+180)
% show at once. Shows |B| contour + in-plane (Br,Bz) quiver on that cut.
launchCrossSectionViewer(R, gridExtentR, gridExtentZ, ID_opening, OD_opening, depth_opening);

%% LOCAL FUNCTIONS
function f = windingCurrent(region, ~, x0, y0, Jmag)
% Azimuthal current density about local axis (x0,y0); sign(Jmag) sets polarity
[TH, ~, ~] = cart2pol(region.x - x0, region.y - y0, region.z);
f = Jmag * [-sin(TH); cos(TH); zeros(size(TH))];
end

function plotColoredStreamlines(verts, X, Y, Z, CGrid)
% surface()+EdgeColor="interp" colors a 3-D line by CGrid; plot3/streamline can't
for i = 1:numel(verts)
    v = verts{i};
    if size(v, 1) < 2
        continue;
    end
    x = v(:, 1)'; y = v(:, 2)'; z = v(:, 3)';
    c = interp3(X, Y, Z, CGrid, x, y, z, 'linear');
    surface([x; x], [y; y], [z; z], [c; c], ...
        FaceColor="none", EdgeColor="interp", LineWidth=1.5);
end
end

function launchCrossSectionViewer(R, Rmax, Zmax, ID_opening, OD_opening, depth_opening)
nR = 50; nZ = 36; % kept modest so each redraw (release-triggered) stays snappy
rg = linspace(-Rmax, Rmax, nR);
zg = linspace(-Zmax, Zmax, nZ);
[Rg, Zg] = meshgrid(rg, zg);

fig = uifigure('Name', 'HET B-field cross-section viewer', 'Position', [100 100 900 650]);
ax = uiaxes(fig, 'Position', [60 130 780 480]);
lbl = uilabel(fig, 'Position', [800 110 100 22], 'Text', '0 / 180 deg');
uilabel(fig, 'Position', [60 60 780 22], 'Text', 'Cross-section angle theta');
sld = uislider(fig, 'Position', [80 80 700 3], 'Limits', [0 180], ...
    'Value', 0, 'MajorTicks', 0:30:180);

channel.ID = ID_opening;
channel.OD = OD_opening;
channel.depth = depth_opening;

updateCrossSectionPlot(ax, R, Rg, Zg, 0, channel);
% Recompute on release (not during drag), since each redraw queries the FEM
% solution directly via interpolateMagneticFlux. lbl flips to "Computing..."
% with a forced drawnow so the delay doesn't look like a stuck UI.
sld.ValueChangedFcn = @(src, ~) onCrossSectionSliderChange(src, ax, R, Rg, Zg, lbl, channel);
end

function onCrossSectionSliderChange(src, ax, R, Rg, Zg, lbl, channel)
lbl.Text = 'Computing...';
drawnow;
thetaDeg = src.Value;
updateCrossSectionPlot(ax, R, Rg, Zg, thetaDeg, channel);
lbl.Text = sprintf('%.0f / %.0f deg', thetaDeg, thetaDeg + 180);
end

function updateCrossSectionPlot(ax, R, Rg, Zg, thetaDeg, channel)
theta = deg2rad(thetaDeg);
X = Rg * cos(theta);
Y = Rg * sin(theta);
Z = Zg;

BI = R.interpolateMagneticFlux(X, Y, Z);
Bx = reshape(BI.Bx, size(Rg));
By = reshape(BI.By, size(Rg));
Bz = reshape(BI.Bz, size(Rg));
Br = Bx*cos(theta) + By*sin(theta); % in-plane radial component of this cut
Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);

cla(ax);
contourf(ax, Rg, Zg, Bmag, 30, LineStyle="none");
hold(ax, 'on');
colormap(ax, 'parula');
cb = colorbar(ax); cb.Label.String = '|B| (T)';

skip = max(1, floor(size(Rg, 1)/20));
quiver(ax, Rg(1:skip:end, 1:skip:end), Zg(1:skip:end, 1:skip:end), ...
    Br(1:skip:end, 1:skip:end), Bz(1:skip:end, 1:skip:end), 'k');

% Channel outline: annulus [ID/2, OD/2] x [-depth/2, depth/2] on both
% half-planes of this cut (channel is assumed centered on z = 0, same as
% the coils)
rectangle(ax, Position=[channel.ID/2, -channel.depth/2, ...
    (channel.OD - channel.ID)/2, channel.depth], ...
    EdgeColor="w", LineStyle="--", LineWidth=1.5);
rectangle(ax, Position=[-channel.OD/2, -channel.depth/2, ...
    (channel.OD - channel.ID)/2, channel.depth], ...
    EdgeColor="w", LineStyle="--", LineWidth=1.5);

xlabel(ax, 'r (m)  [negative r = opposite half-plane]');
ylabel(ax, 'z (m)');
title(ax, sprintf('B field cross-section at theta = %.0f / %.0f deg', thetaDeg, thetaDeg + 180));
axis(ax, 'equal');
hold(ax, 'off');
end
