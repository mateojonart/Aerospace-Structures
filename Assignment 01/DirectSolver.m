classdef DirectSolver < handle

    properties (Access = private)
        k
        p
        data
        f1
        extForces
    end

    methods (Access = public)
        function obj = DirectSolver(cParams)
            obj.k = cParams.K;
            obj.p = cParams.p;
            obj.f1 = cParams.f;
            obj.data = cParams.data;
            obj.extForces = cParams.F;
        end

        function [u,r] = solve(obj)
            [up,vp] = obj.applyBC();
            obj.f1 = obj.pointLoads();
            [u,r] = obj.solveSystem(up,vp);
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

        function f = pointLoads(obj)
            ni = obj.data.ni;
            Fext = zeros(size(obj.f1));
            Fext(nod2dof(ni,obj.extForces(:,1),obj.extForces(:,2))) = obj.extForces(:,3);
            f = obj.f1 + Fext;
        end

        function [u,r] = solveSystem(obj,up,vp)
            ndof = obj.data.ndof;
            vf = setdiff((1:ndof)',vp);
            u = zeros(ndof,1);
            u(vp) = up;
            u(vf) = obj.k(vf,vf)\(obj.f1(vf)-obj.k(vf,vp)*u(vp));
            r = obj.k(vp,:)*u - obj.f1(vp);
        end
    end
end