function [x_path,Tn_path] = BeamDiscretization(b,step)

    n_el = b / step;
    x_path = zeros(1,n_el+1);
    Tn_path = zeros(n_el,2);

    aux = 1;
    for ee = 0:(n_el-1)
        x_path(aux) = step * ee;
        Tn_path(aux,:) = [aux aux+1];
        aux = aux + 1;
    end
    x_path(n_el+1) = b;

end

