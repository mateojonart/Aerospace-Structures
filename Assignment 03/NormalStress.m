function [sigma,s_sigma] = NormalStress(data,x_path,Tn_path,Centroid,Inertias,Moments)
    s_sigma = zeros(2,data.nel);
    sigma = zeros(2,data.nel);

    Ix = Inertias(1) - Inertias(3)^2/Inertias(2);
    Iy = Inertias(2) - Inertias(3)^2/Inertias(1);
    Mx = Moments(1) + Moments(2)*Inertias(3)/Inertias(2);
    My = Moments(2) + Moments(1)*Inertias(3)/Inertias(1);

    for ii = 1:data.nel
        P1 = x_path(Tn_path(ii,1),:);
        P2 = x_path(Tn_path(ii,2),:);
        l = norm(P2-P1);
        if ii > 1
            s_sigma(1,ii) = s_sigma(2,ii-1);
        end
        s_sigma(2,ii) = s_sigma(1,ii) + l;
        sigma(1,ii) = (P1(2) - Centroid(2))*Mx/Ix - (P1(1) - Centroid(1))*My/Iy;
        sigma(2,ii) = (P2(2) - Centroid(2))*Mx/Ix - (P2(1) - Centroid(1))*My/Iy;        
    end
end