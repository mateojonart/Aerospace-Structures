function f = ElementForce(e,x,Tn,fe,me)    
    % Get element's length
    l = abs(x(Tn(e,2)) - x(Tn(e,1)));
    % Contribution of distributed shear load
    f_b = fe * l * [1/2 l/12 0 1/2 -l/12 0]';
    % Contribution of distributed torsion
    f_t = me * l * [0 0 1/2 0 0 1/2]';
    % Compute total element force vector
    f = f_b + f_t;
end

