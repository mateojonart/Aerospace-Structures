function [K,f] = assemblyFunction(data,Td,Kel,Fel)
    K = zeros(data.ndof,data.ndof);
    f = zeros(data.ndof,1);
    for ii = 1:data.nel
        for jj = 1:(data.nne*data.ni)
            f(Td(ii,jj)) = f(Td(ii,jj)) + Fel(jj,ii);
            for kk = 1:(data.nne*data.ni)
                K(Td(ii,jj),Td(ii,kk)) = K(Td(ii,jj),Td(ii,kk)) + Kel(jj,kk,ii);
            end
        end
    end 
end