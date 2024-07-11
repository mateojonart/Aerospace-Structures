function K = ElementStiffness(e,x,Tn,Tm,m)
    
    % Get element's length
    l = abs(x(Tn(e,2)) - x(Tn(e,1)));
    % Get element's material properties
    E = m(Tm(e),1);
    G = m(Tm(e),2);
    I = m(Tm(e),3);
    J = m(Tm(e),4);
    % Compute element's stifness matrix
    K_b = E * I / l^3 * [  12    6*l  0  -12    6*l  0;
                          6*l  4*l^2  0   -6*l  2*l^2  0;
                            0      0  0      0      0  0;
                          -12   -6*l  0     12   -6*l  0;
                          6*l  2*l^2  0   -6*l  4*l^2  0;
                            0      0  0      0      0  0];
    % Compute element's stifness matrix for torsion
    K_t = G * J / l * [0  0   0  0  0   0;
                       0  0   0  0  0   0;
                       0  0   1  0  0  -1;
                       0  0   0  0  0   0;
                       0  0   0  0  0   0;
                       0  0  -1  0  0   1];
    % Compute total element's stifness matrix
    K = K_b + K_t;

end

