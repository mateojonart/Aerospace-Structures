classdef VfComputer < handle

    properties (Access = private)
        data
        prescribDOF
    end

    methods (Access = public)

        function obj = VfComputer(cParams)
            obj.data = cParams.data;
            obj.prescribDOF = cParams.prescribDOF;
        end

        function [up,vp,vf] = compute(obj)
            [up,vp] = obj.applyBC();
            vf = obj.calculateVf(vp);
        end

    end

    methods (Access = private)

        function [up,vp] = applyBC(obj)
            ni = obj.data.ni;
            up = zeros(size(obj.prescribDOF,1),1);
            vp = zeros(size(obj.prescribDOF,1),1);
            for ii = 1:size(obj.prescribDOF,1)
                vp(ii) = nod2dof(ni,obj.prescribDOF(ii,1),obj.prescribDOF(ii,2));
                up(ii) = obj.prescribDOF(ii,3);
            end
        end

        function vf = calculateVf(obj,vp)
            vf = setdiff((1:obj.data.ndof)',vp);
        end

    end

end