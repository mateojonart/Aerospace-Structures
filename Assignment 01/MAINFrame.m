%% STRUCTURAL PROBLEM CODE STRUCTURE

clear
close all
clc

%% 1) PREPROCESS

% 1.1 Input data (define your input parameters here)
data.ni = 2;  % Degrees of freedom per node
young = 71e9; % Young's modulous
grav = 9.81;
geo.d1 = 36e-3;
geo.d2 = 30e-3;
geo.d3 = 20e-3;
geo.t1 = 1.5e-3;
geo.t2 = 1.2e-3;
geo.t3 = 1e-3;

% 1.2 Build geometry (mesh)
% Nodal coordinates matrix
coord = [% column_1 = x-coord , column_2 = y-coord , ...    
        0        0
    0.459   -0.054
    1.125        0
    0.315    0.486
    0.864    0.486
];
data.nnod = size(coord,1); % Number of nodes 
data.nd = size(coord,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom

% Nodal connectivities matrix
connec = [% column_1 = element node 1 , column_2 = element node 2, ...
    1 4
    1 2
    2 4
    4 5
    2 5
    3 5
];
data.nel = size(connec,1); % Number of elements 
data.nne = size(connec,2); % Number of nodes in a bar

% Create degrees of freedom connectivities matrix
connecDOF = connectDOF(data,connec);

% Material properties matrix
matProp = [% Each column corresponds to a material property (area, Young's modulus, etc.)
    young pi*(((geo.d1+geo.t1)/2)^2-((geo.d1-geo.t1)/2)^2) 0 pi/4*(((geo.d1+geo.t1)/2)^4-((geo.d1-geo.t1)/2)^4)
    young pi*(((geo.d2+geo.t2)/2)^2-((geo.d2-geo.t2)/2)^2) 0 pi/4*(((geo.d2+geo.t2)/2)^4-((geo.d2-geo.t2)/2)^4)
    young pi*(((geo.d3+geo.t3)/2)^2-((geo.d3-geo.t3)/2)^2) 0 pi/4*(((geo.d3+geo.t3)/2)^4-((geo.d3-geo.t3)/2)^4)
];

% Material connectivities matrix
connecMat = [% Each row is the material (row number in 'm') associated to each element
    3
    3
    2
    2
    1
    1
];

% 1.3 Input boundary conditions
% Fixed nodes matrix
prescribDOF = [% Each row is a prescribed degree of freedom | column_1 = node, column_2 = direction, column_3 = value of prescribed displacement
    1 1 0
    1 2 0
    3 1 0
    3 2 0
];

% Point loads matrix
extForces = [% Each row is a point force component | column_1 = node, column_2 = direction (1 = x-direction, 2 = y-direction), column_3 = force magnitude
    2 2 -0.45*grav*75
    4 2 -0.5*grav*75
    5 2 -0.05*grav*75
    5 1 75*2.5
];

%% 2) SOLVER

s.data = data;
s.matProp = matProp;
s.connecDOF = connecDOF;
s.connec = connec;
s.connecMat = connecMat;
s.coord = coord;
s.prescribDOF = prescribDOF;
s.extForces = extForces;
s.type = 'iterative';
s.tol = 1e-6;
s.maxIt = 15000;

s.K = GlobalStiffnessMatrixComputer(s).compute();
s.f = GlobalForceVectorComputer(s).compute();
[A,b] = SystemCreator(s).create();
s.x0 = zeros(size(b));
u2 = zeros(size(b,1),1);
[up2,vp2,vf2] = VfComputer(s).compute();
u2(vp2) = up2;
solver = Solver.create(s);
u2(vf2) = solver.solve(A,b);

% 2.1.1 Compute element stiffness matrices
Kel = stiffnessFunction(data,coord,connec,matProp,connecMat);

% 2.1.2 Compute element force vectors
Fel = forceFunction(data,coord,connec,matProp,connecMat); 

% 2.2 Assemble global stiffness matrix
[K,f] = assemblyFunction(data,connecDOF,Kel,Fel);

% 2.3.1 Apply prescribed DOFs
[up,vp] = applyBC(data,prescribDOF);

% 2.3.2 Apply point loads
f = pointLoads(data,f,extForces);

% 2.4 Solve system
[u,r] = solveSystem(data,K,f,up,vp);

% 2.5 Compute stress
sig = stressFunction(data,coord,connec,matProp,connecMat,connecDOF,u);

%% 3) POSTPROCESS

scale = 1; % Set a number to visualize deformed structure properly
units = 'MPa'; % Define in which units you're providing the stress vector
safetyfactor = 2.5; % Set the convenient safety factor

% plot2DBars(data,x,Tn,u,sig*1e-6,scale,units);

% Returns which elements do fail by buckling
Fail = buckling(data,coord,connec,sig,matProp,connecMat);