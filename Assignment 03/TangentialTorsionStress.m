function [tau_torsion,s_tau_torsion] = TangentialTorsionStress(data,x_path,Tn_path,m,Tm,Inertias,TorsionMom,Ain)
    s_tau_torsion = zeros(2,data.nel);
    tau_torsion = zeros(2,data.nel);

    for ii = 1:data.nel
        t = m(Tm(ii),1);
        P1 = x_path(Tn_path(ii,1),:);
        P2 = x_path(Tn_path(ii,2),:);
        l = norm(P2 - P1);
        if ii > 1
            s_tau_torsion(1,ii) = s_tau_torsion(2,ii-1);
        end
        s_tau_torsion(2,ii) = s_tau_torsion(1,ii) + l;
        if data.open
            tau_torsion(:,ii) = TorsionMom * t / Inertias(4);
        else
            tau_torsion(:,ii) = TorsionMom/(2*Ain*t);
        end
    end
end