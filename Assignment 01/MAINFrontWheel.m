%% STRUCTURAL PROBLEM CODE STRUCTURE

clear
close all
clc

%% 1) PREPROCESS

% 1.1 Input data (define your input parameters here)
data.ni = 2;  % Degrees of freedom per node
ERim = 70e9; % Young's modulous [Pa]
ARim = 140/1000^2; % Surface [m^^2]
IRim = 1470/1000^4; % Inertia [mm^^4]
ESpokes = 210e9; % Young's modulous [Pa]
ASpokes = 3.8/1000^2; % Surface [m^^2]
ISpokes = 1.15/1000^4; % Inertia [m^4]

% 1.2 Build geometry (mesh)
% Nodal coordinates matrix [m]
x = 0.35.*[% column_1 = x-coord , column_2 = y-coord , ...    
        0        0
    cosd(22.5)  sind(22.5)
    cosd(22.5)  -sind(22.5)
    sind(22.5)  -cosd(22.5)
    -sind(22.5) -cosd(22.5)
    -cosd(22.5) -sind(22.5)
    -cosd(22.5) sind(22.5)
    -sind(22.5) cosd(22.5)
    sind(22.5)  cosd(22.5)
];
data.nnod = size(x,1); % Number of nodes 
data.nd = size(x,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom

% Nodal connectivities matrix
Tn = [% column_1 = element node 1 , column_2 = element node 2, ...
    9 2
    2 3
    3 4
    4 5
    5 6
    6 7
    7 8
    8 9
    1 2
    1 3
    1 4
    1 5
    1 6
    1 7
    1 8
    1 9
];
data.nel = size(Tn,1); % Number of elements 
data.nne = size(Tn,2); % Number of nodes in a bar

% Create degrees of freedom connectivities matrix
Td = connectDOF(data,Tn);

% Material connectivities matrix
Tm = [% Each row is the material (row number in 'm') associated to each element
    1
    1
    1
    1
    1
    1
    1
    1
    2
    2
    2
    2
    2
    2
    2
    2
];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AQUÍ
% 1.3 Input boundary conditions
% Fixed nodes matrix
p = [% Each row is a prescribed degree of freedom | column_1 = node, column_2 = direction, column_3 = value of prescribed displacement
    4 2 0
    5 1 0
    5 2 0
];

% Point loads matrix [N] --> Reaccions del programa anterior del xassís.
F = [% Each row is a point force component | column_1 = node, column_2 = direction (1 = x-direction, 2 = y-direction), column_3 = force magnitude
    1 1 186.53525
    1 2 -347.3415
    ];

sigma0 = 110616239; % Starting point 
safety = 5; % Initialization (put here any value not close to 2.5)

%while abs(safety - 2.5) >= 1e-6

% Material properties matrix
m = [% Each column corresponds to a material property (area, Young's modulus, etc.)
    ERim ARim 0 IRim 
    ESpokes ASpokes sigma0 ISpokes
];

% 2) SOLVER

% 2.1.1 Compute element stiffness matrices
Kel = stiffnessFunction(data,x,Tn,m,Tm);

% 2.1.2 Compute element force vectors
Fel = forceFunction(data,x,Tn,m,Tm); 

% 2.2 Assemble global stiffness matrix
[K,f] = assemblyFunction(data,Td,Kel,Fel);

% 2.3.1 Apply prescribed DOFs
[up,vp] = applyBC(data,p);

% 2.3.2 Apply point loads
f = pointLoads(data,f,F);

% 2.4 Solve system
[u,r] = solveSystem(data,K,f,up,vp);

% 2.5 Compute stress
sig = stressFunction(data,x,Tn,m,Tm,Td,u);

% 3) POSTPROCESS

scale = 100; % Set a number to visualize deformed structure properly
units = 'MPa'; % Define in which units you're providing the stress vector
safetyfactor = 2.5; % Set the convenient safety factor

% plot2DBars(data,x,Tn,u,sig*1e-6,scale,units);

% Returns which elements do fail by buckling
[SigCrit,Fail] = buckling(data,x,Tn,sig,m,Tm);

safety = min(abs(SigCrit))/abs(min(sig(9:end)));
sigma0 = sigma0 + 1; % Adjust here the step for more or less precision

%end 
