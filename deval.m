function [SxInt, SpxInt] = deval(sol, xInt)
% DEVAL  Evaluate ODE solution and its derivative (Octave compatibility).
%   [SxInt, SpxInt] = deval(sol, xInt) evaluates the solution sol at points
%   xInt. Returns interpolated values SxInt and derivatives SpxInt.
%
%   This is a simplified replacement for MATLAB's deval using cubic
%   interpolation for values and finite differences for derivatives.

    t = sol.x;
    y = sol.y;

    % Interpolate solution values
    SxInt = zeros(size(y,1), length(xInt));
    for i = 1:size(y,1)
        SxInt(i,:) = interp1(t, y(i,:), xInt, 'pchip');
    end

    % Compute derivatives via finite differences on the dense output
    if nargout > 1
        SpxInt = zeros(size(SxInt));
        for i = 1:size(y,1)
            SpxInt(i,:) = gradient(SxInt(i,:), xInt);
        end
    end
end
