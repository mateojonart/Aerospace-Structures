function [u,r] = solveSystem(data,K,f,up,vp)

    % Determine the free DOFs indices vector
    vf = setdiff((1:data.ndof)', vp);
    % Initialize the global DOFs vector
    u = zeros(data.ndof, 1);
    % Assign restrictions to the restricted DOFs on the global DOFs vector
    u(vp) = up;
    % Compute the free DOFs on the global DOFs vector
    tic
    u(vf) = (K(vf,vf)\eye(size(vf,1),size(vf,1))) * (f(vf) - K(vf,vp)*u(vp));
    toc
    % Compute the reactions loads of the restricted DOFs
    r = K(vp,:)*u - f(vp);

end