function [tau_shear,s_tau_shear] = TangentialShearStress(data,x_path,Tn_path,m,Tm,Centroid,Inertias,Shears,ShearCenter,Ain)
    s_tau_shear = zeros(2,data.nel);
    tau_shear = zeros(2,data.nel);

    if data.open
        q = 0;
    else
        q = (Shears(2)*(Centroid(1) - ShearCenter(1)) - Shears(1)*(Centroid(2) - ShearCenter(2))) / (2*Ain);
    end

    Ix = Inertias(1) - Inertias(3)^2/Inertias(2);
    Iy = Inertias(2) - Inertias(3)^2/Inertias(1);
    Sx = Shears(1) - Shears(2)*Inertias(3)/Inertias(1);
    Sy = Shears(2) - Shears(1)*Inertias(3)/Inertias(2);

    for ii = 1:data.nel
        t = m(Tm(ii),1);
        P1 = x_path(Tn_path(ii,1),:);
        P2 = x_path(Tn_path(ii,2),:);
        deltax = P2(1) - P1(1);
        deltay = P2(2) - P1(2);
        l = norm(P2 - P1);
        if ii > 1
            s_tau_shear(1,ii) = s_tau_shear(2,ii-1);
        end
        s_tau_shear(2,ii) = s_tau_shear(1,ii) + l;
        tau_shear(1,ii) = q / t;
        q = q - Sx*t*l*(deltax/2 + P1(1) - Centroid(1)) / Iy - Sy*t*l*(deltay/2 + P1(2) - Centroid(2)) / Ix;
        tau_shear(2,ii) = q / t;
    end
end