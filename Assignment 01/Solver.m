classdef Solver < handle

    properties (Access = public)
        solverSelector
        cParams
    end

    methods (Access = public)
        function obj = Solver(solverSelector,cParams)
            obj.solverSelector = solverSelector;
            obj.cParams = cParams;
        end

        function [u,r] = solve(obj)
            switch obj.solverSelector
                case 1                                      % Direct
                    thisObj = DirectSolver(obj.cParams);
                    [u,r] = thisObj.solve;
                case 2                                      % Iterative
                    thisObj = IterativeSolver(obj.cParams);
                    [u,r] = thisObj.solve;
            end
        end
    end
end