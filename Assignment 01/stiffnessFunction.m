function [Kel,xel] = stiffnessFunction(data,x,Tn,m,Tm)
    Kel = zeros(data.nne*data.ni,data.nne*data.ni,data.nel);
    for ii = 1:data.nel
        [xel] = [x(Tn(ii,:),:)];
        l = (sqrt((xel(2,1) - xel(1,1))^2 + (xel(2,2) - xel(1,2))^2));
        c = (xel(2,1)-xel(1,1))/l;
        s = (xel(2,2)-xel(1,2))/l;
        E = m(Tm(ii),1);
        A = m(Tm(ii),2);    
        Kel(:,:,ii)=E*A/l.*[c^2 c*s -c^2 -c*s; 
                           c*s s^2 -c*s -s^2;
                           -c^2 -c*s c^2 c*s; 
                           -c*s -s^2 c*s s^2];
    end
end