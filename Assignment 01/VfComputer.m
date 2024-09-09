classdef VfComputer < handle

    properties (Access = private)
        data
        p
    end

    methods (Access = public)

        function obj = VfComputer(cParams)
            obj.data = cParams.data;
            obj.p = cParams.prescribDOF;
        end

        function [up,vp,vf] = compute(obj)
            [up,vp] = obj.applyBC();
            vf = obj.calculateVf(vp);
        end

    end

    methods (Access = private)

        function [up,vp] = applyBC(obj)
            ni = obj.data.ni;
            up = zeros(size(obj.p,1),1);
            vp = zeros(size(obj.p,1),1);
            for ii = 1:size(obj.p,1)
                vp(ii) = nod2dof(ni,obj.p(ii,1),obj.p(ii,2));
                up(ii) = obj.p(ii,3);
            end
        end

        function vf = calculateVf(obj,vp)
            vf = setdiff((1:obj.data.ndof)',vp);
        end

    end

end