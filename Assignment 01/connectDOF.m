function Td = connectDOF(data,Tn)
    Td = zeros(data.nel,data.nne*data.ni);
    for ii = 1:data.nel
        aa = 1;
        for jj = 1:data.nne
            for kk = 1:data.ni
                Td(ii,aa) = nod2dof(data.ni,Tn(ii,jj),kk);
                aa = aa + 1;
            end      
        end
    end
end