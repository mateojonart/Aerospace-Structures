classdef SystemCreator < handle

    properties (Access = private)
        K
        prescribDOF
        data
        f
    end

    methods (Access = public)
        function obj = SystemCreator(cParams)
            obj.K = cParams.K;
            obj.prescribDOF = cParams.prescribDOF;
            obj.f = cParams.f;
            obj.data = cParams.data;
        end

        function [A,b] = create(obj)
            [up,vp] = obj.applyBC();
            [A,b] = obj.solveSystem(up,vp);
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

        function [A,b] = solveSystem(obj,up,vp)
            ndof = obj.data.ndof;
            vf = setdiff((1:ndof)',vp);
            u = zeros(ndof,1);
            u(vp) = up;
            A = obj.K(vf,vf);
            b = obj.f(vf)-obj.K(vf,vp)*u(vp);
        end
    end
end