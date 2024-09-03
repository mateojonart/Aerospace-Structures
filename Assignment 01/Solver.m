classdef Solver < handle

    methods (Static)

        function [u,r] = create(cParams)
            switch cParams.type
                case {"direct"}
                    [u,r] = DirectSolver();
                case {"iterative"}
                    [u,r] = IterativeSolver();
            end
        end
    end
end
