clear
close all
clc

%% 1.1 Input data (define your input parameters here)
data.ni = 2;  % Degrees of freedom per node
ERim = 70e9; % Young's modulous [Pa]
ARim = 140/1000^2; % Surface [m^^2]
IRim = 1470/1000^4; % Inertia [mm^^4]
ESpokes = 210e9; % Young's modulous [Pa]
ASpokes = 3.8/1000^2; % Surface [m^^2]
ISpokes = 1.15/1000^4; % Inertia [m^4]

% 1.2 Build geometry (mesh)
% Nodal coordinates matrix [m]
coord = 0.35.*[% column_1 = x-coord , column_2 = y-coord , ...    
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
data.nnod = size(coord,1); % Number of nodes 
data.nd = size(coord,2);   % Problem dimension
data.ndof = data.nnod*data.ni;  % Total number of degrees of freedom

% Nodal connectivities matrix
connec = [% column_1 = element node 1 , column_2 = element node 2, ...
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
data.nel = size(connec,1); % Number of elements 
data.nne = size(connec,2); % Number of nodes in a bar

% Create degrees of freedom connectivities matrix
connecDOF = connectDOF(data,connec);

% Material connectivities matrix
connecMat = [% Each row is the material (row number in 'm') associated to each element
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

% 1.3 Input boundary conditions
% Fixed nodes matrix
prescribDOF = [% Each row is a prescribed degree of freedom | column_1 = node, column_2 = direction, column_3 = value of prescribed displacement
    4 2 0
    5 1 0
    5 2 0
];

% Point loads matrix [N] --> Reaccions del programa anterior del xassís.
extForces = [% Each row is a point force component | column_1 = node, column_2 = direction (1 = x-direction, 2 = y-direction), column_3 = force magnitude
    1 1 0.964750000000009
    1 2 -388.4085
    ];

sigma0 = 45257350; % Starting point (must be below the exact valueyou're tying to find)
safety = 5; % Initialization (put here any value not close to 2.5)

% Material properties matrix
matProp = [% Each column corresponds to a material property (area, Young's modulus, etc.)
    ERim ARim 0 IRim 
    ESpokes ASpokes sigma0 ISpokes
];

% CParams Construction
s.data = data;
s.matProp = matProp;
s.connecDOF = connecDOF;
s.connec = connec;
s.connecMat = connecMat;
s.coord = coord;
s.prescribDOF = prescribDOF;
s.extForces = extForces;
s.type = 'iterative';
s.tol = 1e-8;
s.maxIt = 15000;
s.x0 = zeros(s.data.ndof-size(prescribDOF,1),1);

%% Results to then compare:

s.trueDisplacements = importdata("RearWheelDisplacements.mat");