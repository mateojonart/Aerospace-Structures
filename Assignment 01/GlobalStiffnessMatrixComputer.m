classdef GlobalStiffnessMatrixComputer < handle

    properties (Access = private)
        data
        matProp
        coord
        connec
        connecMat
        connecDOF
    end

    methods (Access = public)

        function obj = GlobalStiffnessMatrixComputer(cParams)
            obj.data = cParams.data;
            obj.matProp = cParams.matProp;
            obj.coord = cParams.coord;
            obj.connec = cParams.connec;
            obj.connecMat = cParams.connecMat;
            obj.connecDOF = cParams.connecDOF;
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
                [xel] = [obj.coord(obj.connec(ii,:),:)];
                l = (sqrt((xel(2,1) - xel(1,1))^2 + (xel(2,2) - xel(1,2))^2));
                c = (xel(2,1)-xel(1,1))/l;
                s = (xel(2,2)-xel(1,2))/l;
                E = obj.matProp(obj.connecMat(ii),1);
                A = obj.matProp(obj.connecMat(ii),2);
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
                        K(obj.connecDOF(ii,jj),obj.connecDOF(ii,kk)) = K(obj.connecDOF(ii,jj),obj.connecDOF(ii,kk)) + Kel(jj,kk,ii);
                    end
                end
            end
        end
    end
end