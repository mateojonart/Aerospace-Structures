function [x_el,S_el,Mb_el,Mt_el] = ElementInternalForces(x,Tn,Td,Kel,u,data)

    % Initialization of vectors
    x_el = zeros(2,data.nel);
    S_el = zeros(2,data.nel);
    Mb_el = zeros(2,data.nel);
    Mt_el = zeros(2,data.nel);

    % Get element's vertices coords.
    x_el(:,:) = x(Tn(:,:))';

    % Loop over each element
    for ee = 1:data.nel
        % Initialization of element's vertices DOFs vector
        u_el = zeros(data.nne*data.ni,1);
        % Loop over each DOF of the element
        for ii = 1:(data.nne*data.ni)
            % Assign the corresponding DOF value
            u_el(ii) = u(Td(ee,ii));
        end
        % Compute internal forces
        fint_el = Kel(:,:,ee) * u_el;
        % Assign cross-section loads
        S_el(1,ee) = -fint_el(1);
        S_el(2,ee) = fint_el(4);
        Mb_el(1,ee) = -fint_el(2);
        Mb_el(2,ee) = fint_el(5);
        Mt_el(1,ee) = -fint_el(3);
        Mt_el(2,ee) = fint_el(6);
    end

end

