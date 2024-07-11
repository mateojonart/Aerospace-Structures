function [Centroid,ShearCenter,Atot,Inertias,Ain] = SectionProperties(data,x_path,Tn_path,m,Tm)

    % Centroid calculation:
    l = zeros(data.nel,1);
    A = zeros(data.nel,1);
    ElemCentroid = zeros(data.nel,2);
    Atot = 0;
    Centroid = [0 0];
    for ii = 1:data.nel
        t = m(Tm(ii),1);
        P1 = x_path(Tn_path(ii,1),:);
        P2 = x_path(Tn_path(ii,2),:);
        l(ii) = norm(P2-P1);
        A(ii) = l(ii)*t;
        ElemCentroid(ii,:) = (P1 + P2)/2;
        Atot = Atot + A(ii);
        Centroid = Centroid + A(ii)*x_path(ii,:);
    end
    Centroid = Centroid/Atot;

    % Inertia calculation:
    Ixx = 0;
    Iyy = 0;
    Ixy = 0;
    J = 0; 
    Ain = 0;
    for jj = 1:data.nel
        t = m(Tm(jj),1);
        P1 = x_path(Tn_path(jj,1),:);
        P2 = x_path(Tn_path(jj,2),:);
        deltax = P2(1) - P1(1);
        deltay = P2(2) - P1(2);
        Ixx = Ixx + A(jj)*deltay^2/12 + A(jj)*(ElemCentroid(jj,2) - Centroid(2))^2;
        Iyy = Iyy + A(jj)*deltax^2/12 + A(jj)*(ElemCentroid(jj,1) - Centroid(1))^2;
        Ixy = Ixy + A(jj)*deltax*deltay/12 + A(jj)*(ElemCentroid(jj,1) - Centroid(1))*(ElemCentroid(jj,2) - Centroid(2));
        if data.open
            J = J + l(jj) * t^3 / 3;
        else
            Ain = Ain + norm(cross([P1(1)-Centroid(1),P1(2)-Centroid(2),0],[deltax,deltay,0]))/2;
            J = J + l(jj)/t;
        end
    end
    
    if data.open == false
        J = 4*Ain^2/J;
    end
    Inertias = [Ixx Iyy Ixy J];

    % Shear center calculation:
    q = [0 0];
    ShearCenter = Centroid;
    for kk = 1:data.nel
        t = m(Tm(kk),1);
        P1 = x_path(Tn_path(kk,1),:);
        P2 = x_path(Tn_path(kk,2),:);
        deltax = P2(1) - P1(1);
        deltay = P2(2) - P1(2); 
        A1 = (Iyy*deltay/2 - Ixy*deltax/2) / (Ixx*Iyy - Ixy^2);
        A2 = (Ixx*deltax/2 - Ixy*deltay/2) / (Ixx*Iyy - Ixy^2); 
        B1 = (Iyy*(P1(2)-Centroid(2)) - Ixy*(P1(1)-Centroid(1))) / (Ixx*Iyy - Ixy^2);
        B2 = (Ixx*(P1(1)-Centroid(1)) - Ixy*(P1(2)-Centroid(2))) / (Ixx*Iyy - Ixy^2);
        C = (P1(1)-Centroid(1))*deltay - (P1(2)-Centroid(2))*deltax;
        ShearCenter = ShearCenter + C*(q - t*l(kk)*([A1 -A2]/3 + [B1 -B2]/2));
        q = q - t*l(kk)*([A1 -A2] + [B1 -B2]);
    end
end 
