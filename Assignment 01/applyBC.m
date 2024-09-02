function [up,vp] = applyBC(data,p)
    up = zeros(size(p,1),1);
    vp = zeros(size(p,1),1);
    for ii = 1:size(p,1)
        vp(ii) = nod2dof(data.ni,p(ii,1),p(ii,2));
        up(ii) = p(ii,3);
    end
end