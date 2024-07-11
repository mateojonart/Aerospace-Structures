function [up,vp] = applyBC(data,p)

    vp = zeros(size(p,1), 1);  % Restricted DOF indexes
    up = zeros(size(p,1), 1);  % Corresponding DOF restriction

    % Loop each row of matrix p
    for ii = 1:size(p,1)
        % Determine the DOF index
        vp(ii) = nod2dof(data.ni,p(ii,1),p(ii,2));
        % Determine the corresponding DOF restriction
        up(ii) = p(ii,3);
    end

end