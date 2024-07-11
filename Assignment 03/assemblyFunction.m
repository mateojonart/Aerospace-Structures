% Function to assemble the global stiffness matrix and the force element
% vector.

function [K,f] = assemblyFunction(data,Td,Kel,fel)

    K = zeros(data.ndof, data.ndof);
    f = zeros(data.ndof, 1);

    % Loop each element
    for ee = 1:data.nel
        % Loop each DOF of each element
        for ii = 1:(data.nne*data.ni)
            % Add the element force vector to the global force vector
            f(Td(ee,ii)) = f(Td(ee,ii)) + fel(ii,ee);
            % Loop each DOF of each element
            for jj = 1:(data.nne*data.ni)
                % Add the element stiffness matrix to the gloabl stiffness
                % matrix
                K(Td(ee,ii),Td(ee,jj)) = K(Td(ee,ii),Td(ee,jj)) + Kel(ii,jj,ee);
            end
        end
    end

end