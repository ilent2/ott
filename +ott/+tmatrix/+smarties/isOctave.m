function [ bIsOctave ] = isOctave()
%% isOctave
% Test if we are running octave
%
%   Returns true if we're running octave, false otherwise
%
% Dependency: 
% none

    bIsOctave = exist('OCTAVE_VERSION', 'builtin')~=0;
end

