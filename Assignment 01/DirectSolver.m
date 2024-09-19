classdef DirectSolver < Solver

    methods (Static, Access = protected)

        function x = solve(A,b)
            x = A\b;
        end

    end

    methods (Access = public)
        function obj = DirectSolver(cParams)
            obj@Solver(cParams);
        end
    end

end