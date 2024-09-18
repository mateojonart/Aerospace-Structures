classdef TestDisp < handle

    properties (Access = public)
        K
        f
        data
        matProp
        connecDOF
        connec
        connecMat
        coord
        prescribDOF
        extForces
        type
        tol
        maxIt
        x0
    end

    properties (Access = private)
        trueDisplacements
    end

    methods (Access = public)

        function obj = TestDisp(cParams)
            obj.init(cParams);
        end

        function  u = solveAndCheckProblem(obj,sParams)
            u = obj.solveProblem(sParams);
            obj.verifyResults(u)
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.data = cParams.data;
            obj.matProp = cParams.matProp;
            obj.connecDOF = cParams.connecDOF;
            obj.connec = cParams.connec;
            obj.coord = cParams.coord;
            obj.prescribDOF = cParams.prescribDOF;
            obj.extForces = cParams.extForces;
            obj.type = cParams.type;
            obj.tol = cParams.tol;
            obj.maxIt = cParams.maxIt;
            obj.x0 = cParams.x0;
            obj.connecMat = cParams.connecMat;
            obj.trueDisplacements = cParams.trueDisplacements;
        end

        function u = solveProblem(obj,sParams)
            obj.K = GlobalStiffnessMatrixComputer(sParams).compute();
            obj.f = GlobalForceVectorComputer(sParams).compute();
            [A,b] = SystemCreator(obj).create();
            u = zeros(obj.data.ndof,1);
            [up,vp,vf] = VfComputer(sParams).compute();
            u(vp) = up;
            solver = Solver.create(sParams);
            u(vf) = solver.solve(A,b);
        end

        function verifyResults(obj,u)
            err = norm(obj.trueDisplacements-u);
            if err <= 1e-12
                disp('Test passed.');
            else
                disp(['Test failed. The error is ', num2str(err),'.']);
            end
        end

    end
end