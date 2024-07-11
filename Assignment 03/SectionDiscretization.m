function [x_path,Tn_path,Tm_path] = SectionDiscretization(x,Tn,path_step)

    aux = 1;
    for ii = 1:size(Tn,1)
        P1 = x(Tn(ii,1),:);
        P2 = x(Tn(ii,2),:);
        l = norm(P2-P1);
        n = floor(l/path_step);
        real_step = l / n;
        direction_vec = (P2 - P1)/l;

        for jj = 0:n-1
            x_path(aux,:) = P1 + direction_vec * jj*real_step;
            Tn_path(aux,:) = [aux aux+1];
            Tm_path(aux,1) = ii;
            aux = aux + 1;
            disp(jj)
        end
        disp(ii);
    end
    x_path(end+1,:) = x(1,:);
end