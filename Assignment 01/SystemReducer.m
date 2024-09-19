classdef SystemReducer < handle

    properties (Access = private)
        data
        uprescrib
        vprescrib
        vfree
    end

    methods (Access = public)
        function obj = SystemReducer(cParams)
            obj.init(cParams);
        end

        function [A,b] = reduceSystem(obj,K,f)
            ndof = obj.data.ndof;
            vf   = obj.vfree;
            vp   = obj.vprescrib;
            up   = obj.uprescrib;
            u = zeros(ndof,1);
            u(vp) = up;
            A = K(vf,vf);
            b = f(vf)-K(vf,vp)*u(vp);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.uprescrib   = cParams.up;
            obj.vprescrib   = cParams.vp;
            obj.vfree       = cParams.vf;
            obj.data        = cParams.data;
        end
    end
end