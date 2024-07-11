clc
clear
close all

%% SECTION ANALYSIS

% INPUT data_section

% Open/Closed section
data_section.open = 0;

% Geometry of the cross section:
geo.c = 2; % [m]
geo.xip = 0.3*geo.c;
geo.d = 0.3*geo.c;
geo.h1 = 0.25*geo.c;
geo.h2 = 0.15*geo.c;
geo.t1 = 22e-3;
geo.t2 = 15e-3;
geo.t3 = 3.5e-3;

% Nodal coordinates matrix before discretization:
x = [        geo.d                       geo.h2/2;
             geo.d     geo.h2+((geo.h1-geo.h2)/2);
                 0                         geo.h2;
                 0                              0;
             geo.d             -(geo.h1-geo.h2)/2];

% Nodal connectivities matrix before discretization:
Tn = [1 2;
      2 3;
      3 4;
      4 5;
      5 1];

% Material properties matrix:
m_section = [geo.t1;
             geo.t3;
             geo.t2;
             geo.t3
             geo.t1];

% DISCRETIZATION OF THE CROSS SECTION

% Distance between points of the path:
path_step = 0.00005; % [m]

% Process of discretization in a uniform path:
[x_path,Tn_path,Tm_path] = SectionDiscretization(x,Tn,path_step);

% Show resulting figure:
% scatter(x_path(:,1),x_path(:,2),'Marker','o')

% Problem section data
data_section.ni = 1; % Degrees of freedom per node
data_section.nnod = size(x_path,1); % Number of nodes
data_section.nel = length(Tm_path); % Number of elements
data_section.ndof = data_section.ni*data_section.nnod; % Total number of DOFs

% COMPUTATION OF OF THE SECTION PROPERTIES
[Centroid,ShearCenter,Atot,Inertias,Ain] = SectionProperties(data_section,x_path,Tn_path,m_section,Tm_path);
Inertias(3) = 0;

% Is closed section, shear center is considered to be in the centroid
if data_section.open ~= true
    ShearCenter = Centroid;
end

%% SECTION UNITARY LOADS ANALYSIS

% Bending moments vector: [-Mz My]
BendingMoms = [-1 0];

% Shears vector: [-Tz Ty]
Shears = [0 1];

% Torsion moment: [Mx]
TorsionMom = -1;

% Stresses calculation:
[sigma,s_sigma] = NormalStress(data_section,x_path,Tn_path,Centroid,Inertias,BendingMoms);
[tau_shear,s_tau_shear] = TangentialShearStress(data_section,x_path,Tn_path,m_section,Tm_path,Centroid,Inertias,Shears,ShearCenter,Ain);
[tau_torsion,s_tau_torsion] = TangentialTorsionStress(data_section,x_path,Tn_path,m_section,Tm_path,Inertias,TorsionMom,Ain);

% Show the results:
% PlotSection(x_path,Tn_path,sigma/1e3,'KPa');
% PlotSection(x_path,Tn_path,tau_shear/1e3,'KPa');
% PlotSection(x_path,Tn_path,tau_torsion/1e3,'KPa'); 

%% BEAM ANALYSIS

% Geometry of the wing:
geo.b = 16; % [m]
geo.za = 0.25*geo.c; % [m]
geo.zm = 0.48*geo.c; % [m]
geo.be = 0.25*geo.b; % [m]
geo.ze = 0.3*geo.c; % [m]
geo.me = 2.1e3; % [kg]
geo.lambda = 140; % [kg/m]

% Physical properties:
geo.E = 210e9; % [Pa]
geo.G = 80e9; % [Pa]
geo.g = 9.8065; % [N/kg] or [m/s^2]
geo.rho_inf = 1.225; % [kg/m^3]
geo.v_inf = 750/3.6; % [m/s]
geo.cl = 0.1;

% Material properties matrix: [E G Izz J]
m_beam = [geo.E geo.G Inertias(1) Inertias(4)];

% Number of panels along the beam:
% NumPanels = [1 2 4 8 16 32 64 128 256 512 1024 2048]; % To study the convergence
NumPanels = 512; % Reference number of panels

% Space allocation for vectors
ElementVertForce = cell(length(NumPanels),1); % Element equivalent vertical force
ElementEquivalentTorsionMoment = cell(length(NumPanels),1); % Element equivalent torsion moment
ElementCenterCoord = cell(length(NumPanels),1); % Element center coordinates
VerticalDeflection = cell(length(NumPanels),1); % Nodes' vertical deflection
BendingRotation = cell(length(NumPanels),1); % Nodes' bending angle
TwistRotation = cell(length(NumPanels),1); % Nodes' twist angle
x_beam_plot = cell(length(NumPanels),1); % Nodes' coordinates
Shear_beam = cell(length(NumPanels),1); % Internal element shear forces
BendingMoment_beam = cell(length(NumPanels),1); % Internal element bending moment 
TorsionalMoment_Beam = cell(length(NumPanels),1); % Internal element torsional moment
TipDeflectionError = zeros(length(NumPanels),1); % Tip Deflection relative error
TipTwistError = zeros(length(NumPanels),1); % Tip twist relative error

% Loop over every number of panels
for pp = 1:length(NumPanels)

    % Beam discretization step
    step = geo.b / NumPanels(pp); % Always returns a round number since we choose appropriate numbers of panels 
    
    % Beam discretization
    [x_beam,Tn_beam] = BeamDiscretization(geo.b,step);
    x_beam_plot{pp} = x_beam;
    
    % Problem beam data
    data_beam.ni = 3; % (vertical deflection , bending rotation , torsional rotation)
    data_beam.nel = size(Tn_beam,1); % Number of elements
    data_beam.nne = size(Tn_beam,2); % Number of nodes per element
    data_beam.nnod = length(x_beam); % Number of nodes 
    data_beam.ndof = data_beam.nnod * data_beam.ni; % Total number of DOFs
    
    % Connectivities matrix for the beam
    Td_beam = connectDOF(data_beam,Tn_beam);

    % Material properties vector of the beam
    Tm_beam = ones(data_beam.nel,1);
    
    % Prescribed DOFs
    p = [1 1 0; % Deflection at node 1 = 0 
         1 2 0; % Bending at node 1 = 0
         1 3 0;]; % Torsion at node 1 = 0
    
    % First panel where the engine weight is applied (second one will be n+1)
    n = round(geo.be / step);

    % Point loads matrix
    F = [n+1 1                                               -geo.me*geo.g;
         n+1 3          -geo.me*geo.g*(geo.xip+geo.d-ShearCenter(1)-geo.ze)];
    
    % BEAM LOADS ANALYSIS
    
    % Space allocation for vectors and matrices
    Kel = zeros(data_beam.nne*data_beam.ni,data_beam.nne*data_beam.ni,data_beam.nel); % Stiffness matrices of each element
    fel = zeros(data_beam.nne*data_beam.ni,data_beam.nel); % Equivalent force on each element
    
    % Lift distribution on each node
    l = 0.5 * geo.rho_inf * geo.v_inf^2 * geo.c * geo.cl .* sqrt(1 - (x_beam(:) ./ geo.b).^2);

    % Space allocation for vectors and matrices used for plotting later:
    feplot = zeros(data_beam.nel,1); % Vertical force on each element
    meplot = zeros(data_beam.nel,1); % Equivalent torsion moment on each element
    xplot = zeros(data_beam.nel,1); % Coordinates of the center of each element
    
    % Loop over each element to study its loads
    for ee = 1:data_beam.nel
    
        % Compute element stiffness matrices
        Kel(:,:,ee) = ElementStiffness(ee,x_beam,Tn_beam,Tm_beam,m_beam);
    
        % Compute element force vectors
        fe_l = (l(Tn_beam(ee,1)) + l(Tn_beam(ee,2))) / 2; % Vertical force due to lift (positive)
        fe_w = -geo.lambda * geo.g; % Vertical force due to weight (negative)
    
        % Resulting force:
        fe = fe_l + fe_w;
        feplot(ee) = fe;
    
        % Equivalent torsion moment about the shear center
        me = fe_l * (geo.xip + geo.d - ShearCenter(1) - geo.za) + fe_w * (geo.xip + geo.d - ShearCenter(1) - geo.zm); 
        meplot(ee) = me;

        % Calculation of the center of the element
        xplot(ee) = step * (ee - 1) + step / 2;

        % Computation of force element vector
        fel(:,ee) = ElementForce(ee,x_beam,Tn_beam,fe,me);

    end

    % Element center coordinates:
    ElementCenterCoord{pp} = xplot;

    % External forces and moments over each element:
    ElementVertForce{pp} = feplot; % Vertical force
    ElementEquivalentTorsionMoment{pp} = meplot; % Torsion moment
    
    % Assembly of the global stiffness matrix
    [K,f] = assemblyFunction(data_beam,Td_beam,Kel,fel);

    % Apply the point load of the engine
    f = pointLoads(data_beam,f,F);

    % Apply restricted DOFs
    [up,vp] = applyBC(data_beam,p);

    % Solve the system of equations inverting the free DOFs matrix 
    [u,~] = solveSystem(data_beam,K,f,up,vp);

    % Compute internal forces of each DOF
    [~,S_el,Mb_el,Mt_el] = ElementInternalForces(x_beam,Tn_beam,Td_beam,Kel,u,data_beam);
    Shear_beam{pp} = [S_el(1,:) S_el(2,end)]; % Internal shear
    BendingMoment_beam{pp} = [Mb_el(1,:) Mb_el(2,end)]; % Internal bending moment
    TorsionalMoment_Beam{pp} = [Mt_el(1,:) Mt_el(2,end)]; % Internal torsion moment

    % Vertical deflection:
    VerticalDeflection{pp} = u(1:3:data_beam.ndof); % DOFs 1,4,7,10,...

    % Bending Rotation:
    BendingRotation{pp} = u(2:3:data_beam.ndof); % DOFs 2,5,8,11,...

    % Torsional Rotation:
    TwistRotation{pp} = u(3:3:data_beam.ndof); % DOFs 3,6,9,12,...

end

%% CONVERGENCE ANALYSIS:

% Loop over each number of panels to compute errors
for ii = 1:length(NumPanels)
    % Tip deflection error
    TipDeflectionError(ii) = abs((VerticalDeflection{ii}(end) - VerticalDeflection{end}(end))./ VerticalDeflection{end}(end)) * 100;
    % Tip twist error
    TipTwistError(ii) = abs((TwistRotation{ii}(end) - TwistRotation{end}(end))./ TwistRotation{end}(end)) * 100;
end

%% MOST CRITICAL POINT ANALYSIS (VOM MISES CRITERION):

% Space allocation for equivalent tensions and its indexes
sigmaVM_section = zeros(length(x_beam),1); % Maximum equivalent tension of each section
idx_sigmaVM_section = zeros(length(x_beam),1); % Its index in the x_path vector

% Loop over each section of the beam (each node)
for ii = 1:length(x_beam)

    % Computation of normal stress
    [sigma,s_sigma] = NormalStress(data_section,x_path,Tn_path,Centroid,Inertias,[-BendingMoment_beam{1}(ii) 0]);
    sigma_vec = [sigma(1,:) sigma(2,end)];

    % Computation of tangential shear stress
    [tau_shear,s_tau_shear] = TangentialShearStress(data_section,x_path,Tn_path,m_section,Tm_path,Centroid,Inertias,[0 Shear_beam{1}(ii)],ShearCenter,Ain);
    tau_shear_vec = [tau_shear(1,:) tau_shear(2,end)];

    % Computation of tangential torsion stress
    [tau_torsion,s_tau_torsion] = TangentialTorsionStress(data_section,x_path,Tn_path,m_section,Tm_path,Inertias,TorsionalMoment_Beam{1}(ii),Ain);
    tau_torsion_vec = [tau_torsion(1,:) tau_torsion(2,end)];

    % Total tangential stress (in modulous)
    tau = abs(tau_shear_vec) + abs(tau_torsion_vec);

    % Finding the maximum value and its index inside x_path vector
    [sigmaVM_section(ii), idx_sigmaVM_section(ii)] = max(sqrt(sigma_vec.^2 + 3*tau.^2));

    if ii == 129
        sigmaVM_section_crit = sqrt(sigma_vec.^2 + 3*tau.^2);
    end

    disp(ii);

end

% Find the maximum of the maximums and its index within the x_beam array:
[sigmaVM_beam,idx_sigmaVM_beam] = max(sigmaVM_section./1e6); % In MPa

% Position of the critical point
posMax_beam = x_beam(idx_sigmaVM_beam);
posMax_section = x_path(idx_sigmaVM_section(idx_sigmaVM_beam),:);

%% PLOTS

% figure; hold on; grid on;
% for ii = 1:length(NumPanels)
%     plot(ElementCenterCoord{ii}, ElementVertForce{ii},'LineWidth',1.5);
% end
% title('Vertical external force along the wing','FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('x (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Vertical Force (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% % legend('4 Panels','8 Panels','16 Panels','32 Panels','64 Panels','128 Panels','256 Panels','512 Panels','1024 Panels','2048 Panels','Location','northwest','NumColumns',2,'FontSize',11.5);
% ax = gca;
% ax.LineWidth = 1.5;
% hold off;
% 
% figure; hold on; grid on;
% for ii = 1:length(NumPanels)
%     plot(ElementCenterCoord{ii}, ElementEquivalentTorsionMoment{ii},'LineWidth',1.5);
% end
% title('Equivalent torsion along the wing','FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('x (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Equivalent torsion (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% % legend('4 Panels','8 Panels','16 Panels','32 Panels','64 Panels','128 Panels','256 Panels','512 Panels','1024 Panels','2048 Panels','Location','northwest','NumColumns',2,'FontSize',11.5);
% ax = gca;
% ax.LineWidth = 1.5;
% hold off;

% figure; hold on; grid on;
% for ii = 1:length(NumPanels)
%     plot(x_beam_plot{ii}, VerticalDeflection{ii},'LineWidth',1.5);
% end
% title('Vertical Deflection along the wing','FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('Wing position (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Vertical Deflection u_y (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% % legend('4 Panels','8 Panels','16 Panels','32 Panels','64 Panels','128 Panels','256 Panels','512 Panels','1024 Panels','2048 Panels','Location','northwest','NumColumns',2,'FontSize',11.5);
% ax = gca;
% ax.LineWidth = 1.5;
% hold off;
% 
% figure; hold on; grid on;
% for ii = 1:length(NumPanels)
%     plot(x_beam_plot{ii}, BendingRotation{ii},'LineWidth',1.5);
% end
% title('Section bending rotation along the wing','FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('x (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Section rotation angle (rad)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% % legend('4 Panels','8 Panels','16 Panels','32 Panels','64 Panels','128 Panels','256 Panels','512 Panels','1024 Panels','2048 Panels','Location','southeast','NumColumns',2,'FontSize',11.5);
% ax = gca;
% ax.LineWidth = 1.5;
% hold off;

% figure; hold on; grid on;
% for ii = 1:length(NumPanels)
%     plot(x_beam_plot{ii}, TwistRotation{ii},'LineWidth',1.5);
% end
% title('Section twist angle along the wing','FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('x (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Section twist angle (rad)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% % legend('4 Panels','8 Panels','16 Panels','32 Panels','64 Panels','128 Panels','256 Panels','512 Panels','1024 Panels','2048 Panels','Location','southeast','NumColumns',2,'FontSize',11.5);
% ax = gca;
% ax.LineWidth = 1.5;
% hold off;
% 
% figure; hold on; grid on;
% for ii = 1:length(NumPanels)
%     plot(x_beam_plot{ii}, Shear_beam{ii},'LineWidth',1.5);
% end
% title('Shear force along the wing','FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('x (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Shear Force (N)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% % legend('4 Panels','8 Panels','16 Panels','32 Panels','64 Panels','128 Panels','256 Panels','512 Panels','1024 Panels','2048 Panels','Location','northwest','NumColumns',2,'FontSize',11.5);
% ax = gca;
% ax.LineWidth = 1.5;
% hold off;

% figure; hold on; grid on;
% for ii = 1:length(NumPanels)
%     plot(x_beam_plot{ii}, BendingMoment_beam{ii},'LineWidth',1.5);
% end
% title('Bending moment along the wing','FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('x (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Bending moment (N·m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% % legend('4 Panels','8 Panels','16 Panels','32 Panels','64 Panels','128 Panels','256 Panels','512 Panels','1024 Panels','2048 Panels','Location','northwest','NumColumns',2,'FontSize',11.5);
% ax = gca;
% ax.LineWidth = 1.5;
% hold off;
% 
% figure; hold on; grid on;
% for ii = 1:length(NumPanels)
%     plot(x_beam_plot{ii}, TorsionalMoment_Beam{ii},'LineWidth',1.5);
% end
% title('Torsional moment along the wing','FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('x (m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Torsional moment (N·m)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% % legend('4 Panels','8 Panels','16 Panels','32 Panels','64 Panels','128 Panels','256 Panels','512 Panels','1024 Panels','2048 Panels','Location','northwest','NumColumns',2,'FontSize',11.5);
% ax = gca;
% ax.LineWidth = 1.5;
% hold off;

% figure(1);
% subplot(2,1,1);
% semilogx(NumPanels,TipTwistError,'LineWidth',1.5);
% grid on;
% title("Relative error of wing tip's twist",'FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('Number of panels','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Relative error (%)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% xlim([1 max(NumPanels)]);
% set(gca,'LineWidth',1.5); 
% 
% subplot(2,1,2);
% semilogx(NumPanels,TipDeflectionError,'LineWidth',1.5);
% grid on;
% title("Relative error of wing tip's deflection",'FontSize',16,'FontWeight','bold','FontAngle','italic');
% xlabel('Number of panels','FontSize',14,'FontWeight','bold','FontAngle','italic');
% ylabel('Relative error (%)','FontSize',14,'FontWeight','bold','FontAngle','italic');
% xlim([1 max(NumPanels)]);
% set(gca,'LineWidth',1.5);
