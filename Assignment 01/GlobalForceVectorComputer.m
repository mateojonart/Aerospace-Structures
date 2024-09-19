classdef GlobalForceVectorComputer < handle

    properties (Access = private)
        data
        matProp
        coord
        connec
        connecMat
        connecDOF
        extForces
    end

    methods (Access = public)

        function obj = GlobalForceVectorComputer(cParams)
            obj.data = cParams.data;
            obj.matProp = cParams.matProp;
            obj.coord = cParams.coord;
            obj.connec = cParams.connec;
            obj.connecMat = cParams.connecMat;
            obj.connecDOF = cParams.connecDOF;
            obj.extForces = cParams.extForces;
        end

        function f = compute(obj)
            Fel = obj.calculateElemForce();
            f = obj.assembleElemForces(Fel);
        end
    end

    methods (Access = private)

        function Fel = calculateElemForce(obj)
            nne = obj.data.nne;
            ni = obj.data.ni;
            nel = obj.data.nel;
            Fel = zeros(nne*ni,nel);
            for ii = 1:nel
                xel = [obj.coord(obj.connec(ii,:),:)];
                l = sqrt((xel(2,1) - xel(1,1))^2 + (xel(2,2) - xel(1,2))^2);
                c = (xel(2,1)-xel(1,1))/l;
                s = (xel(2,2)-xel(1,2))/l;
                R = [c s 0 0;
                    -s c 0 0;
                    0 0 c s;
                    0 0 -s c];
                A = obj.matProp(obj.connecMat(ii),2);
                sigma0 = obj.matProp(obj.connecMat(ii),3);
                Fel(:,ii) = -sigma0.*A.*R'*[-1; 0; 1; 0];
            end
        end

        function f = assembleElemForces(obj,Fel)
            ndof = obj.data.ndof;
            nel = obj.data.nel;
            nne = obj.data.nne;
            ni = obj.data.ni;
            f_prov = zeros(ndof,1);
            for ii = 1:nel
                for jj = 1:(nne*ni)
                    f_prov(obj.connecDOF(ii,jj)) = f_prov(obj.connecDOF(ii,jj)) + Fel(jj,ii);
                end
            end
            Fext = zeros(size(f_prov));
            Fext(nod2dof(ni,obj.extForces(:,1),obj.extForces(:,2))) = obj.extForces(:,3);
            f = f_prov + Fext;
        end

    end
end