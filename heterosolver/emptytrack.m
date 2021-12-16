function [ trjempty ] = emptytrack()
%EMPTYTRACK generate an empty track, no molecule or background intensity
%   Usage: [trjempty] = emptytrack()

trjempty.no = zeros(0,1);
trjempty.track = zeros(0,5);

end

