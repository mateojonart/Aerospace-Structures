function Fel = forceFunction(data,x,Tn,m,Tm)
    Fel = zeros(data.nne*data.ni,data.nel);
    for ii = 1:data.nel
        [xel] = [x(Tn(ii,:),:)];
        l = sqrt((xel(2,1) - xel(1,1))^2 + (xel(2,2) - xel(1,2))^2);
        c = (xel(2,1)-xel(1,1))/l;
        s = (xel(2,2)-xel(1,2))/l;
        R = [c s 0 0; 
             -s c 0 0;
             0 0 c s; 
             0 0 -s c];
        A = m(Tm(ii),2);
        sigma0 = m(Tm(ii),3);
        Fel(:,ii)=-sigma0.*A.*R'*[-1; 0; 1; 0];
    end
end