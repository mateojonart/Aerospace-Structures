classdef Solver < handle

    methods (Static)

        function stype = create(cParams)
            switch cParams.type
                case {"direct"}
                    stype = Direct_Solver();
                case {"iterative"}
                    stype = Iterative_pcg_Solver();
                otherwise 
                    error('Invalid solver type')
            end
        end
    end
end
