classdef TestDisp < handle

    properties (Access = private)
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
        trueDisplacements
    end

    methods (Access = public)

        function obj = TestDisp(cParams)
            obj.init(cParams);
        end

        function  u = solveAndCheckProblem(obj)
            u = obj.solveProblem();
            obj.checkResults(u)
        end
    end

    methods (Access = private)

        function u = solveProblem(obj)
            [K,f] = obj.computeStiffnessAndForce();
            solver = obj.createSolver();
            u = solver.solveSystem(K,f);
        end

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

        function [K,f] = computeStiffnessAndForce(obj)
            sKC.data      = obj.data;
            sKC.coord     = obj.coord;
            sKC.connec    = obj.connec;
            sKC.matProp   = obj.matProp;
            sKC.connecMat = obj.connecMat;
            sKC.connecDOF = obj.connecDOF;
            sKC.extForces = obj.extForces;
            stiffnessComputer = GlobalStiffnessMatrixComputer(sKC);
            K = stiffnessComputer.compute();
            forceComputer = GlobalForceVectorComputer(sKC);
            f = forceComputer.compute();
        end

        function solver = createSolver(obj)
            sS.x0          = obj.x0;
            sS.tol         = obj.tol;
            sS.data        = obj.data;
            sS.type        = obj.type;
            sS.maxIt       = obj.maxIt;
            sS.prescribDOF = obj.prescribDOF;
            solver = Solver.create(sS);
        end

        function checkResults(obj,u)
            err = norm(obj.trueDisplacements-u);
            if err <= 1e-12
                disp('Test passed.');
            else
                disp(['Test failed. The error is ', num2str(err),'.']);
            end
        end
    end
end