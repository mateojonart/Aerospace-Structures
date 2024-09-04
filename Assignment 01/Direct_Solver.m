classdef Direct_Solver < handle

    methods (Static)

        function x = solve(A,b)
            x = A\b;
        end

    end

end