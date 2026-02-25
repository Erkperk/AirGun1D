function [dq, dBubble] = airgunRHS(q, t, bubble, deps)
% AIRGUNRHS  Right-hand side for the coupled airgun-bubble system.
%   Standalone function for Octave compatibility (no nested functions in classdef).

    flowState = deps.flowStateR(q);

    if t >= deps.physConst.AirgunCutoffTime || flowState == deps.SUBSONIC_INFLOW
        dq = q.*0;
        dBubble = bubbleRHS(bubble, 0, 0, 0, 0, 0, deps.physConst);
        return
    end

    p_b = bubblePressure(bubble, deps.physConst);
    dq = deps.D(q) + deps.closure_l(q);

    if flowState == deps.SUBSONIC_OUTFLOW
        dq = dq + deps.closure_r_out_sub(q, p_b);
    elseif flowState == deps.SUPERSONIC_OUTFLOW
        % No bc required
    else
        q_r = deps.e_R' * q;
        fprintf('q_r(1) => %g\n', q_r(1));
        fprintf('q_r(2) => %g\n', q_r(2));
        fprintf('q_r(3) => %g\n', q_r(3));
        error('Undefined behaviour. flowState = %d, t = %f', flowState, t)
    end

    qHat = deps.e_R' * q;
    rho_a = qHat(1);
    v_a = qHat(2) / qHat(1);
    e_a = qHat(3);
    p_a = deps.p_fun(qHat);

    % turbulent energy dissipation
    C = 0;
    gun_area = 12.5 * 6.4516e-4; % [m^2] (12.5 in^2)
    port_area = 12 * 6.4516e-4;  % [m^2] (12 in^2)
    deltaP = C * rho_a * v_a^2 * (gun_area/port_area)^2;
    dBubble = bubbleRHS(bubble, rho_a, v_a, e_a, p_a - deltaP, deps.A, deps.physConst);
end
