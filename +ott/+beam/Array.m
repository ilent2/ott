classdef Array < ott.beam.ArrayType
% A array of generic beams.
% Inherits from :class:`ArrayType`.
%
% This class stores a internal list of :class:`Beam` instances.
% Data from coherent beams is combined after the field calculation
% methods are called.  Data for incoherent beams is combined by
% the :class:`ArrayType` class methods.
%
% For now, the array content should all have the same optical
% frequency and medium refractive index.  May change in future.
%
% ``omega`` and ``index_medium`` are updated when the data is changed.
% Changing these properties does not change the beam content.
%
% Properties
%   - data        -- Array of coherent/incoherent beams

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    data      % Array of coherent/incoherent beams
  end

  properties (Hidden, SetAccess=protected)
    dataInternal
  end

  methods
    function beam = Array(varargin)
      % Construct a new beam array instance.
      %
      % Usage
      %   beam = Array(data, ...)
      %
      % Parameters
      %   - data (ott.beam.Beam) -- Array of beams.
      %
      % Unmatched parameters are passed to base class.

      p = inputParser;
      p.addOptional('data', ott.beam.Beam.empty(), ...
          @(x) isa(x, 'ott.beam.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = ott.beam.ArrayType(unmatched{:});
      beam.data = p.Results.data;
    end

    function [E, vswfData] = efield(beam, xyz, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, vswfData] = beam.efield(xyz, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y;z].
      %
      % See :class:`Beam` for additional information.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData());
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      vswfData = p.Results.data;

      E = zeros([size(xyz), numel(beam.data)]);
      for ii = 1:numel(beam.data)
        [E(:, :, ii), vswfData] = beam.data(ii).efield(xyz, ...
            unmatched{:}, 'data', vswfData);
      end

      if strcmpi(beam.arrayType, 'coherent')
        E = sum(E, 3);
      end
    end

    function [H, vswfData] = hfield(beam, xyz, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Usage
      %   [H, vswfData] = beam.hfield(xyz, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y;z].
      %
      % See :class:`Beam` for additional information.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData());
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      vswfData = p.Results.data;

      H = zeros([size(xyz), numel(beam.data)]);
      for ii = 1:numel(beam.data)
        [H(:, :, ii), vswfData] = beam.data(ii).hfield(xyz, ...
            unmatched{:}, 'data', vswfData);
      end

      if strcmpi(beam.arrayType, 'coherent')
        H = sum(H, 3);
      end
    end

    function [E, H, vswfData] = ehfield(beam, xyz, varargin)
      % Calculate the electric and magnetic field (in SI units)
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(xyz, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y;z].
      %
      % See :class:`Beam` for additional information.

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData());
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      vswfData = p.Results.data;

      E = zeros([size(xyz), numel(beam.data)]);
      H = E;
      for ii = 1:numel(beam.data)
        [E(:, :, ii), H(:, :, ii), vswfData] = beam.data(ii).ehfield(xyz, ...
            unmatched{:}, 'data', vswfData);
      end

      if strcmpi(beam.arrayType, 'coherent')
        E = sum(E, 3);
        H = sum(H, 3);
      end
    end

    % TODO: Add force/torque/spin functions to Beam or elsewhere?
    % TODO: Finish other field functions here
  end

  methods % Getters/setters
    function beam = set.data(beam, val)
      assert(isa(val, 'ott.beam.Beam'), ...
          'data must be a array of ott.beam.Beam');

      % Check for incoherent parts
      if strcmpi(beam.arrayType, 'coherent')
        for ii = 1:numel(val)
          if isa(val(ii), 'ott.beam.ArrayType')
            assert(~val(ii).containsIncoherent(), ...
              'Coherent arrays cannot contain incoherent beams');
          end
        end
      end

      % Update index_medium and omega
      try
        beam.index_medium = unique([val.index_medium]);
        beam.omega = unique([val.omega]);
      catch
        error('All sub-beam index-medium/omega must match');
      end

      % Store data
      beam.dataInternal = val;
    end
  end
end
