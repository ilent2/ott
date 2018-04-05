classdef Tmatrix < matlab.mixin.Copyable
%Tmatrix abstract class representing T-matrix of a scattering particle
%
% Tmatrix properties:
%   data          The T-matrix this class encapculates
%
% Tmatrix methods:
%
% Static methods:
%
% Abstract methods:
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

 properties (SetAccess=protected)
  data            % The matrix this class encapsulates
 end

  methods (Static)
    function k_medium = parser_k_medium(p)
      %PARSER_K_MEDIUM helper to get k_medium from a parser object

      if ~isempty(p.Results.k_medium)
        k_medium = p.Results.k_medium;
      elseif ~isempty(p.Results.wavelength_medium)
        k_medium = 2.0*pi/p.Results.wavelength_medium;
      elseif ~isempty(p.Results.index_medium)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_medium');
        end
        k_medium = p.Results.index_medium*2.0*pi/p.Results.wavelength0;
      else
        error('Unable to determine k_medium from inputs');
      end
    end

    function k_particle = parser_k_particle(p)
      %PARSER_K_PARTICLE helper to get k_particle from a parser object

      if ~isempty(p.Results.k_particle)
        k_particle = p.Results.k_particle;
      elseif ~isempty(p.Results.wavelength_particle)
        k_particle = 2.0*pi/p.Results.wavelength_particle;
      elseif ~isempty(p.Results.index_particle)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_particle');
        end
        k_particle = p.Results.index_particle ...
            * 2.0*pi/p.Results.wavelength0;
      else
        error('Unable to determine k_particle from inputs');
      end
    end
  end

 methods (Abstract)
 end

 methods (Access=protected)
  function tmatrix = Tmatrix(data)
    %TMATRIX construct a new T-matrix object
    %
    % TMATRIX() leaves the data uninitialised.
    %
    % TMATRIX(data) initializes the data with the matrix data.

    if nargin == 1
      tmatrix.data = data;
    end
  end
 end

 methods
 end

end
