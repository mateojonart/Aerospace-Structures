classdef GlobalForceVectorComputer < handle

    properties (Access = private)
        data
        m
        x
        tn
        tm
        td
        extForces
    end

    methods (Access = public)

        function obj = GlobalForceVectorComputer(cParams)
            obj.data = cParams.data;
            obj.m = cParams.matProp;
            obj.x = cParams.coord;
            obj.tn = cParams.connec;
            obj.tm = cParams.connecMat;
            obj.td = cParams.connecDOF;
            obj.extForces = cParams.extForces;
        end

        function f = compute(obj)
            Fel = obj.forceFunction();
            f_prov = obj.assemblyFunction(Fel);
            f = obj.pointLoads(f_prov);
        end
    end

    methods (Access = private)

        function Fel = forceFunction(obj)
            nne = obj.data.nne;
            ni = obj.data.ni;
            nel = obj.data.nel;
            Fel = zeros(nne*ni,nel);
            for ii = 1:nel
                xel = [obj.x(obj.tn(ii,:),:)];
                l = sqrt((xel(2,1) - xel(1,1))^2 + (xel(2,2) - xel(1,2))^2);
                c = (xel(2,1)-xel(1,1))/l;
                s = (xel(2,2)-xel(1,2))/l;
                R = [c s 0 0;
                    -s c 0 0;
                    0 0 c s;
                    0 0 -s c];
                A = obj.m(obj.tm(ii),2);
                sigma0 = obj.m(obj.tm(ii),3);
                Fel(:,ii) = -sigma0.*A.*R'*[-1; 0; 1; 0];
            end
        end

        function f = assemblyFunction(obj,Fel)
            ndof = obj.data.ndof;
            nel = obj.data.nel;
            nne = obj.data.nne;
            ni = obj.data.ni;
            f = zeros(ndof,1);
            for ii = 1:nel
                for jj = 1:(nne*ni)
                    f(obj.td(ii,jj)) = f(obj.td(ii,jj)) + Fel(jj,ii);
                end
            end
        end

        function f = pointLoads(obj,f_prov)
            ni = obj.data.ni;
            Fext = zeros(size(f_prov));
            Fext(nod2dof(ni,obj.extForces(:,1),obj.extForces(:,2))) = obj.extForces(:,3);
            f = f_prov + Fext;
        end

    end
end