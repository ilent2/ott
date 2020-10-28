classdef DdaHighMem < ott.tmatrix.dda.Dda
% Slightly faster but more memory intensive version of DDA.
% Inherits from :class:`ott.tmatrix.dda.Dda`.
%
% Properties
%   - interaction     -- Stored interaction matrix
%
% Methods
%   - update_interaction_matrix -- Re-calculate the interaction matrix
%   - interaction_matrix -- Get a usable interaction matrix
%
% Static methods
%   - FromShape       -- Construct from a geometric shape
%
% For other methods/properties, see :class:`Dda`.

  properties (SetAccess=protected)
    interaction     % Stored interaction matrix
  end

  methods (Static)
    function dda = FromShape(shape, varargin)
      % Construct a DDA instance from a geometric shape.
      %
      % Usage
      %   dda = ott.tmatrix.dda.Dda.FromShape(shape, ...)
      %
      % Optional named arguments
      %   - spacing -- (numeric) -- Dipole spacing in wavelength units.
      %     Default: ``1/20``.
      %
      %   - polarizability -- (function_handle | 3x3 numeric) Method to
      %     calculate polarizability or 3x3 tensor for homogeneous material.
      %     Default: ``@(xyz, spacing, ri) polarizability.LDR(spacing, ri)``
      %
      %   - relative_index -- (function_handle | numeric) Method to calculate
      %     relative refractive index or homogeneous value.  Ignored if
      %     polarizability is a 3x3 tensor.
      %
      % For further details and options, see :meth:`Dda.FromShape`.

      dda = ott.tmatrix.dda.Dda.FromShape(shape, varargin{:});
      dda = ott.tmatrix.dda.DdaHighMem(dda);
    end
  end

  methods
    function dda = DdaHighMem(varargin)
      % Construct DDA instance and pre-compute data for interaction matrix.
      %
      % Usage
      %   dda = DdaHigMem(dda)
      %   Convert an existing DDA instance into a pre-computed instance.
      %
      %   dda = DdaHighMem(locations, interaction, ...)
      %   Construct a new DDA instance.  See base class for parameters.

      if nargin == 1
        oldDda = varargin{1};
        assert(isa(oldDda, 'ott.tmatrix.dda.Dda'), ...
            'With only one argument, input must be a Dda instance');
        data = {oldDda.locations, oldDda.polarizability, ...
            'xySymmetry', oldDda.xySymmetry, ...
            'zRotSymmetry', oldDda.zRotSymmetry};
      else
        data = varargin;
      end

      % Construct base
      dda = dda@ott.tmatrix.dda.Dda(data{:});

      % Pre-calculate data
      dda = dda.update_interaction_matrix();
    end

    function dda = update_interaction_matrix(dda)
      % Update the interaction matrix data.
      %
      % This method is called when the class is constructed.
      % If you change the class properties, call this method again.
      %
      % Usage
      %   dda = dda.update_interaction_matrix()

      ott.utils.nargoutCheck(dda, nargout);

      dda.interaction = ott.tmatrix.dda.Dipole.build_field_matrix(...
          dda.locations, dda.locations, ...
          @dda.interaction_efield_column, ...
          'low_memory', false, ...
          'xySymmetry', dda.xySymmetry, ...
          'zRotSymmetry', dda.zRotSymmetry);
    end

    function A = interaction_matrix(dda, varargin)
      % Calculate the interaction matrix from the pre-computed data.
      %
      % This method uses a lot of memory but should be much faster for
      % repeated interaction matrix calculation with rotation/mirror symmetry.
      %
      % Usage
      %   A = dda.intersection_matrix(rorder, parity)
      %
      % Parameters
      %   - parity (enum) -- Parity of incident beam (even or odd).
      %     Only used when using xySymmetry.  Default: ``'even'``.
      %
      %   - rorder (numeric) -- Rotational order of incident beam.
      %     Only used when using zRotSymmetry.  Default: ``0``.

      p = inputParser;
      p.addOptional('parity', 'even');
      p.addOptional('rorder', 0);
      p.parse(varargin{:});

      A = dda.interaction;
      if ~ismatrix(A)
        A = ott.tmatrix.dda.Dipole.combine_rotsym_matrix(...
            A, p.Results.rorder, p.Results.parity);
      end
    end
  end
end

