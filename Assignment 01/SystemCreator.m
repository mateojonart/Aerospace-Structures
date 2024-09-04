classdef SystemCreator < handle

    properties (Access = private)
        k
        p
        data
        f
    end

    methods (Access = public)
        function obj = SystemCreator(cParams)
            obj.k = cParams.K;
            obj.p = cParams.p;
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
            up = zeros(size(obj.p,1),1);
            vp = zeros(size(obj.p,1),1);
            for ii = 1:size(obj.p,1)
                vp(ii) = nod2dof(ni,obj.p(ii,1),obj.p(ii,2));
                up(ii) = obj.p(ii,3);
            end
        end

        function [A,b] = solveSystem(obj,up,vp)
            ndof = obj.data.ndof;
            vf = setdiff((1:ndof)',vp);
            u = zeros(ndof,1);
            u(vp) = up;
            A = obj.k(vf,vf);
            b = obj.f(vf)-obj.k(vf,vp)*u(vp);
        end
    end
end