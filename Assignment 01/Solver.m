classdef Solver < handle

    methods (Static)

        function stype = create(cParams)
            switch cParams.type
                case {"direct"}
                    stype = Direct_Solver();
                case {"iterative"}
                    stype = Iterative_pcg_Solver(cParams);
                otherwise 
                    error('Solver:InvalidType', 'Invalid solver type');
            end
        end
    end
end
