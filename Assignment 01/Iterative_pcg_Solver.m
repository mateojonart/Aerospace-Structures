classdef Iterative_pcg_Solver < handle
    
    properties (Access = private)
        x0
        tol
        maxIt
    end

    methods (Access = public)

        function obj = Iterative_pcg_Solver(cParams)
            Iterative_pcg_Solver.validateParams(cParams);
            obj.x0 = cParams.x0;
            obj.tol = cParams.tol;
            obj.maxIt = cParams.maxIt;
        end

        function x = solve(obj,A,b)
            x = pcg(A,b,obj.tol,obj.maxIt,[],[],obj.x0);
        end

    end

    methods (Static, Access = private)

        function validateParams(cParams)

            if ~isfield(cParams, 'x0') || isempty(cParams.x0)
                error('Iterative_pcg_Solver:MissingField','x0 is required');
            end

            if ~isfield(cParams, 'tol') || isempty(cParams.tol)
                error('Iterative_pcg_Solver:MissingField','tol is required');
            end

            if ~isfield(cParams, 'maxIt') || isempty(cParams.maxIt)
                error('Iterative_pcg_Solver:MissingField','maxIt is required');
            end

            if cParams.tol <= 0
                error('Iterative_pcg_Solver:InvalidTolerance','Tolerance (tol) must be a positive number');
            end

            if ~isnumeric(cParams.maxIt) || numel(cParams.maxIt) ~= 1 || cParams.maxIt <= 0 || cParams.maxIt ~= floor(cParams.maxIt)
                error('Iterative_pcg_Solver:InvalidMaxIterations','Maximum iterations (maxIt) must be a positive integer');
            end

        end

    end

end