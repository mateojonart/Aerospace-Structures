classdef GlobalStiffnessMatrixComputer < handle

    properties (Access = private)
        data
        m
        x
        tn
        tm
        td
    end

    methods (Access = public)

        function obj = GlobalStiffnessMatrixComputer(cParams)
            obj.data = cParams.data;
            obj.m = cParams.m;
            obj.x = cParams.x;
            obj.tn = cParams.tn;
            obj.tm = cParams.tm;
            obj.td = cParams.td;
        end

        function K = compute(obj)
            Kel = obj.stiffnessFunction();
            K = obj.assemblyFunction(Kel);
        end
    end

    methods (Access = private)

        function Kel = stiffnessFunction(obj)
            nne = obj.data.nne;
            ni = obj.data.ni;
            nel = obj.data.nel;
            Kel = zeros(nne*ni,nne*ni,nel);
            for ii = 1:obj.data.nel
                [xel] = [obj.x(obj.tn(ii,:),:)];
                l = (sqrt((xel(2,1) - xel(1,1))^2 + (xel(2,2) - xel(1,2))^2));
                c = (xel(2,1)-xel(1,1))/l;
                s = (xel(2,2)-xel(1,2))/l;
                E = obj.m(obj.tm(ii),1);
                A = obj.m(obj.tm(ii),2);
                Kel(:,:,ii) = E*A/l.*[c^2 c*s -c^2 -c*s;
                    c*s s^2 -c*s -s^2;
                    -c^2 -c*s c^2 c*s;
                    -c*s -s^2 c*s s^2];
            end
        end

        function K = assemblyFunction(obj,Kel)
            ndof = obj.data.ndof;
            nel = obj.data.nel;
            nne = obj.data.nne;
            ni = obj.data.ni;
            K = zeros(ndof,ndof);
            for ii = 1:nel
                for jj = 1:(nne*ni)
                    for kk = 1:(nne*ni)
                        K(obj.td(ii,jj),obj.td(ii,kk)) = K(obj.td(ii,jj),obj.td(ii,kk)) + Kel(jj,kk,ii);
                    end
                end
            end
        end
    end
end