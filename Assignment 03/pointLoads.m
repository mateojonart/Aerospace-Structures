function f = pointLoads(data,f,F)
    Fext = zeros(size(f));
    Fext(nod2dof(data.ni,F(:,1),F(:,2))) = F(:,3);
    f = f + Fext;
end