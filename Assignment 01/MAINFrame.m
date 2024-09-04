%% STRUCTURAL PROBLEM CODE STRUCTURE

clear
close all
clc

%% 1) PREPROCESS

% 1.1 Input data (define your input parameters here)
data.ni = 2;  % Degrees of freedom per node
E = 71e9; % Young's modulous
g = 9.81;
geo.d1 = 36e-3;
geo.d2 = 30e-3;
geo.d3 = 20e-3;
geo.t1 = 1.5e-3;
geo.t2 = 1.2e-3;
geo.t3 = 1e-3;

% 1.2 Build geometry (mesh)
% Nodal coordinates matrix
x = [% column_1 = x-coord , column_2 = y-coord , ...    
        0        0
    0.459   -0.054
    1.125        0
    0.315    0.486
    0.864    0.486
];
data.nnod = size(x,1); % Number of nodes 
data.nd = size(x,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom

% Nodal connectivities matrix
Tn = [% column_1 = element node 1 , column_2 = element node 2, ...
    1 4
    1 2
    2 4
    4 5
    2 5
    3 5
];
data.nel = size(Tn,1); % Number of elements 
data.nne = size(Tn,2); % Number of nodes in a bar

% Create degrees of freedom connectivities matrix
Td = connectDOF(data,Tn);

% Material properties matrix
m = [% Each column corresponds to a material property (area, Young's modulus, etc.)
    E pi*(((geo.d1+geo.t1)/2)^2-((geo.d1-geo.t1)/2)^2) 0 pi/4*(((geo.d1+geo.t1)/2)^4-((geo.d1-geo.t1)/2)^4)
    E pi*(((geo.d2+geo.t2)/2)^2-((geo.d2-geo.t2)/2)^2) 0 pi/4*(((geo.d2+geo.t2)/2)^4-((geo.d2-geo.t2)/2)^4)
    E pi*(((geo.d3+geo.t3)/2)^2-((geo.d3-geo.t3)/2)^2) 0 pi/4*(((geo.d3+geo.t3)/2)^4-((geo.d3-geo.t3)/2)^4)
];

% Material connectivities matrix
Tm = [% Each row is the material (row number in 'm') associated to each element
    3
    3
    2
    2
    1
    1
];

% 1.3 Input boundary conditions
% Fixed nodes matrix
p = [% Each row is a prescribed degree of freedom | column_1 = node, column_2 = direction, column_3 = value of prescribed displacement
    1 1 0
    1 2 0
    3 1 0
    3 2 0
];

% Point loads matrix
F = [% Each row is a point force component | column_1 = node, column_2 = direction (1 = x-direction, 2 = y-direction), column_3 = force magnitude
    2 2 -0.45*g*75
    4 2 -0.5*g*75
    5 2 -0.05*g*75
    5 1 75*2.5
];

%% 2) SOLVER

s.data = data;
s.m = m;
s.td = Td;
s.tn = Tn;
s.tm = Tm;
s.x = x;
s.p = p;
computer = GlobalStiffnessMatrixComputer(s);
[K2,f2] = computer.compute;

% 2.1.1 Compute element stiffness matrices
Kel = stiffnessFunction(data,x,Tn,m,Tm);

% 2.1.2 Compute element force vectors
Fel = forceFunction(data,x,Tn,m,Tm); 

% 2.2 Assemble global stiffness matrix
[K,f] = assemblyFunction(data,Td,Kel,Fel);
s.K = K;
s.f = f;
s.F = F;

% 2.3.1 Apply prescribed DOFs
[up,vp] = applyBC(data,p);

% 2.3.2 Apply point loads
f = pointLoads(data,f,F);
f2 = pointLoads(data,f2,F);

% 2.4 Solve system
[u,r] = solveSystem(data,K,f,up,vp);

% 2.5 Compute stress
sig = stressFunction(data,x,Tn,m,Tm,Td,u);

%% 3) POSTPROCESS

scale = 1; % Set a number to visualize deformed structure properly
units = 'MPa'; % Define in which units you're providing the stress vector
safetyfactor = 2.5; % Set the convenient safety factor

plot2DBars(data,x,Tn,u,sig*1e-6,scale,units);

% Returns which elements do fail by buckling
Fail = buckling(data,x,Tn,sig,m,Tm);