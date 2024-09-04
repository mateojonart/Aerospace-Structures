classdef Solver < handle

    methods (Static)

        function stype = create(cParams)
            switch cParams.type
                case {"direct"}
                    stype = DirectSolver();
                case {"iterative"}
                    stype = IterativeSolver();
                otherwise 
                    error('Invalid solver type')
            end
        end
    end
end
