classdef Solver < handle

    properties (Access = private)
        data
        prescribDOF
    end

    methods (Static)

        function stype = create(cParams)
            switch cParams.type
                case {"direct"}
                    stype = DirectSolver(cParams);
                case {"iterative"}
                    stype = IterativeSolver(cParams);
                otherwise
                    error('Solver:InvalidType', 'Invalid solver type');
            end
        end
    end

    methods (Access = public)

        function obj = Solver(cParams)
            obj.init(cParams);
        end

        function u = solveSystem(obj,K,f)
            a.K           = K;
            a.f           = f;
            a.data        = obj.data;
            a.prescribDOF = obj.prescribDOF;
            computer = DOFComputer(a);
            [up,vp,vf] = computer.computeDOF();
            aS.up   = up;
            aS.vp   = vp;
            aS.vf   = vf;
            aS.data = obj.data;
            system = SystemReducer(aS);
            [KRed,fRed] = system.reduceSystem(K,f);
            u = zeros(obj.data.ndof,1);
            u(vp) = up;
            u(vf) = obj.solve(KRed,fRed);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.data        = cParams.data;
            obj.prescribDOF = cParams.prescribDOF;
        end
    end

    methods (Abstract, Access = protected)
        solve(KRed,fRed);
    end

end
