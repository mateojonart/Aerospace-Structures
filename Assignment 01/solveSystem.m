function [u,r] = solveSystem(data,K,f,up,vp)

    vf = setdiff((1:data.ndof)',vp); % PAS 1
    u = zeros(data.ndof,1); % PAS 2
    u(vp) = up; % PAS 3
    u(vf) = K(vf,vf)\(f(vf)-K(vf,vp)*u(vp)); % PAS 4
    r = K(vp,:)*u - f(vp); % PAS 5
end