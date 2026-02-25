classdef Discretization < handle
    properties
        name         %Short description
        description  %Longer description
        order        %Order of accuracy
    end

    methods
        function printInfo(obj); error('Not implemented'); end
        function n = size(obj); error('Not implemented'); end
        function ts = getTimestepper(obj, opt); error('Not implemented'); end
        function k = getTimestep(obj, opt); error('Not implemented'); end
        function repr = getTimeSnapshot(obj, ts); error('Not implemented'); end
        function [update,hand] = setupPlot(obj, type); error('Not implemented'); end
    end
end