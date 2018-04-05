classdef Bsc < matlab.mixin.Copyable
%Bsc abstract class representing beam shape coefficients
%
% Bsc properties:
%   a               Beam shape coefficients a vector
%   b               Beam shape coefficients b vector
%
% Bsc methods:
%   translateZ      Translates the beam along the z axis
%   translateXyz    Translation to xyz using rotations and z translations
%   translateRtp    Translation to rtp using rotations and z translations
%   farfield        Calculate fields in farfield
%   emfieldXyz      Calculate fields at specified locations
%   set.Nmax        Resize the beam shape coefficient vectors
%   get.Nmax        Get the current size of the beam shape coefficient vectors
%
% Static methods:
%   make_beam_vector    Convert output of bsc_* functions to beam coefficients
%
% Abstract methods:
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
    a           % Beam shape coefficients a vector
    b           % Beam shape coefficients b vector
  end

  methods (Abstract)
  end

  methods (Static)
    % TODO: make_beam_vector
  end

  methods (Access=protected)
    function tmatrix = Bsc()
      %BSC construct a new beam object
    end
  end

  methods

    % TODO: Add a warning when the beam is translated outside nmax2ka(Nmax) 
    % TODO: translateZ
    % TODO: translateXyz
    % TODO: translateRtp
    % TODO: farfield
    % TODO: emfieldXyz
    % TODO: set.Nmax
    % TODO: get.Nmax

  end

end
