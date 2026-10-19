// Dumbbell nanopore geometry: two cylindrical reservoirs joined by a cylindrical pore.
// The membrane is NOT part of the domain; only the electrolyte is.
// Axis of symmetry is z. The membrane occupies |z| < lpn/2.
// All lengths are normalized (pore radius = rpn, normally 1).

SetFactory("OpenCASCADE");

// ---- Parameters (can be overridden with: gmsh pore.geo -setnumber lpn 2) ----
DefineConstant[ rpn = {1,   Name "Pore radius (normalized)"} ];
DefineConstant[ lpn = {1,   Name "Pore length (normalized)"} ];
DefineConstant[ lrn = {5,   Name "Reservoir length (normalized)"} ];
DefineConstant[ rrn = {10,  Name "Reservoir radius (normalized)"} ];
DefineConstant[ ldn = {0.5, Name "Debye length (normalized)"} ];

// ---- Primitives ----
// Cylinder(tag) = {x, y, z,  dx, dy, dz,  radius};  (base center, axis vector, radius)
zb = -lpn/2;   // bottom face of membrane
zt =  lpn/2;   // top face of membrane

Cylinder(1) = {0, 0, zb - lrn,  0, 0, lrn,  rrn};   // bottom reservoir
Cylinder(2) = {0, 0, zt,        0, 0, lrn,  rrn};   // top reservoir
Cylinder(3) = {0, 0, zb,        0, 0, lpn,  rpn};   // pore

// ---- Fuse into one volume (touching faces merge, shared faces disappear) ----
v() = BooleanUnion{ Volume{1}; Delete; }{ Volume{2, 3}; Delete; };

// ---- Physical Groups ----
// Tags change after booleans, so select surfaces geometrically instead of by ID.
eps = 1e-6;
Rmax = rrn + eps;
Ztot = lpn/2 + lrn;

allS[] = Surface{:};                 // every surface in the model
bot[] = {}; top[] = {}; memb_bot[] = {}; memb_top[] = {};
pore_wall[] = {}; wall_bot[] = {}; wall_top[] = {}; other[] = {};

For i In {0:#allS[]-1}
  s = allS[i];
  c[] = CenterOfMass Surface{s};
  z = c[2];
  If (Fabs(z + Ztot) < eps)
    bot[] += s;
  ElseIf (Fabs(z - Ztot) < eps)
    top[] += s;
  ElseIf (Fabs(z - zb) < eps)
    memb_bot[] += s;
  ElseIf (Fabs(z - zt) < eps)
    memb_top[] += s;
  ElseIf (z > zb && z < zt)
    pore_wall[] += s;
  ElseIf (z > -Ztot && z < zb)
    wall_bot[] += s;
  ElseIf (z > zt && z < Ztot)
    wall_top[] += s;
  Else
    other[] += s;
  EndIf
EndFor

Physical Volume("electrolyte") = {v()};
Physical Surface("electrode_bottom") = {bot[]};
Physical Surface("electrode_top")    = {top[]};
Physical Surface("membrane_bottom")  = {memb_bot[]};
Physical Surface("membrane_top")     = {memb_top[]};
Physical Surface("pore_wall")        = {pore_wall[]};
Physical Surface("reservoir_wall_bottom") = {wall_bot[]};
Physical Surface("reservoir_wall_top")    = {wall_top[]};

// ---- Meshing: target element sizes ----
DefineConstant[ fpn  = {0.2, Name "Pore size, fraction of pore radius"} ];
DefineConstant[ frn  = {0.2, Name "Reservoir size, fraction of reservoir radius"} ];
DefineConstant[ rampperc = {0.2, Name "Distance over which size grows outside the pore, fraction of reservoir length"} ];

hdeb = 0.3;         // mesh size in the debeye layer as a fraction of the debeye layer
hpore = fpn * rpn;  // mesh size in the pore
hres  = frn * rrn;  // mesh size in the reservoir

pad   = 1e-3 * rpn;     // extends the box slightly past the pore edges
ramp = rampperc * lrn;

Field[1] = Cylinder;
Field[1].Radius = rpn + pad;
Field[1].VIn  = hpore; Field[1].VOut = hres;
Field[1].ZAxis = lpn + 2*pad;

Field[2] = Min;
Field[2].FieldsList = {1};
Background Field = 2;

Field[3] = Distance;
Field[3].SurfacesList = {pore_wall[], memb_bot[], memb_top[]};
Field[3].Sampling = 100;

ndl = 1;
dramp = 3;

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = hdeb * ldn;          
Field[4].SizeMax = hres;
Field[4].DistMin = ndl * ldn;          
Field[4].DistMax = (ndl + dramp) * ldn;  // ramp out to the reservoir size
Field[4].Sigmoid = 2;

Field[2].FieldsList = {1, 4};          // Min of the pore field and the Debye field

Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;