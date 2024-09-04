classdef Iterative_pcg_Solver < handle
    
    properties (Access = private)
        x0
        tol
        maxIt
    end

    methods (Access = public)

        function obj = Iterative_pcg_Solver(cParams)
            obj.x0 = cParams.x0;
            obj.tol = cParams.tol;
            obj.maxIt = cParams.maxIt;
        end

        function x = solve(obj,A,b)
            x = pcg(A,b,obj.tol,obj.maxIt,[],[],obj.x0);
        end

    end

end