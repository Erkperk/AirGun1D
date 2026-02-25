% Start with all matrix returns. When that works see how we should generalize to non-matrix stuff/nonlinear
classdef Scheme < handle
    properties
        m % Number of points in each direction, possibly a vector
        order % Order accuracy for the approximation
        D % non-stabalized scheme operator
        H % Discrete norm
    end

    methods
        function m = boundary_condition(obj,boundary,type,data); error('Not implemented'); end
        function m = interface(obj,boundary,neighbour_scheme,neighbour_boundary); error('Not implemented'); end
        function N = size(obj); error('Not implemented'); end
    end

    methods(Static)
        % Calculates the matrcis need for the inteface coupling between boundary bound_u of scheme schm_u
        % and bound_v of scheme schm_v.
        %   [uu, uv, vv, vu] = inteface_couplong(A,'r',B,'l')
        function [uu, uv, vv, vu] = interface_coupling(schm_u,bound_u,schm_v,bound_v)
            [uu,uv] = schm_u.interface(bound_u,schm_v,bound_v);
            [vv,vu] = schm_v.interface(bound_v,schm_u,bound_u);
        end
    end
end