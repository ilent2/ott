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

 methods (Abstract)
 end

 methods (Access=protected)
  function tmatrix = Tmatrix()
    %TMATRIX construct a new T-matrix object
  end
 end

 methods
 end

end
