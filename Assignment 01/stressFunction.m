function sig = stressFunction(data,x,Tn,m,Tm,Td,u)
    A = [-1 0 1 0];
    sig = zeros(data.nel, 1);
    for ii = 1:data.nel
        xel = x(Tn(ii,:),:);
        l = sqrt((xel(2,1) - xel(1,1))^2 + (xel(2,2) - xel(1,2))^2);
        c = (xel(2,1)-xel(1,1))/l;
        s = (xel(2,2)-xel(1,2))/l;
        E = m(Tm(ii),1);
        sigma0 = m(Tm(ii),3);
        R = [c s 0 0; -s c 0 0; 0 0 c s; 0 0 -s c];
        uel = u(Td(ii,:),:);
        epsilon = 1/l * A * R * uel;
        sig(ii) = E * epsilon + sigma0;
    end
end
