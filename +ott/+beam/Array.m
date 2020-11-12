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
%   - arrayType   -- Type of array (either coherent/incoherent/array)

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    data      % Array of coherent/incoherent beams
    index_medium    % Medium refractive index
    omega           % Optical frequency
  end

  properties (Hidden, SetAccess=protected)
    dataInternal
  end

  methods
    function bm = Array(varargin)
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

      bm = bm@ott.beam.ArrayType(unmatched{:});
      bm.data = p.Results.data;
    end

    function varargout = efield(beam, varargin)
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

      [varargout{1:nargout}] = beam.singleOutputHelper(...
          @(odata, varargin) odata.efield(varargin{:}), varargin{:});
    end

    function varargout = hfield(beam, varargin)
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

      [varargout{1:nargout}] = beam.singleOutputHelper(...
          @(odata, varargin) odata.hfield(varargin{:}), varargin{:});
    end

    function varargout = ehfield(beam, varargin)
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

      [varargout{1:nargout}] = beam.doubleOutputHelper(...
          @(odata, varargin) odata.ehfield(varargin{:}), varargin{:});
    end

    function varargout = efieldRtp(beam, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, vswfData] = beam.efieldRtp(rtp, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [r;t;p].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.singleOutputHelper(...
          @(odata, varargin) odata.efieldRtp(varargin{:}), varargin{:});
    end

    function varargout = hfieldRtp(beam, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Usage
      %   [H, vswfData] = beam.hfieldRtp(rtp, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [r;t;p].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.singleOutputHelper(...
          @(odata, varargin) odata.hfieldRtp(varargin{:}), varargin{:});
    end

    function varargout = ehfieldRtp(beam, varargin)
      % Calculate the electric and magnetic field (in SI units)
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfieldRtp(rtp, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [r;t;p].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.doubleOutputHelper(...
          @(odata, varargin) odata.ehfieldRtp(varargin{:}), varargin{:});
    end

    function varargout = efarfield(beam, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, vswfData] = beam.efarfield(rtp, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [r;t;p].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.singleOutputHelper(...
          @(odata, varargin) odata.efarfield(varargin{:}), varargin{:});
    end

    function varargout = hfarfield(beam, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Usage
      %   [H, vswfData] = beam.hfarfield(xyz, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [r;t;p].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.singleOutputHelper(...
          @(odata, varargin) odata.hfarfield(varargin{:}), varargin{:});
    end

    function varargout = ehfarfield(beam, varargin)
      % Calculate the electric and magnetic field (in SI units)
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfarfield(rtp, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - rtp (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [r;t;p].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.doubleOutputHelper(...
          @(odata, varargin) odata.ehfarfield(varargin{:}), varargin{:});
    end

    function varargout = eparaxial(beam, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, vswfData] = beam.eparaxial(xy, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - xy (2xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.singleOutputHelper(...
          @(odata, varargin) odata.eparaxial(varargin{:}), varargin{:});
    end

    function varargout = hparaxial(beam, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Usage
      %   [H, vswfData] = beam.hparaxial(xy, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - xy (2xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.singleOutputHelper(...
          @(odata, varargin) odata.hparaxial(varargin{:}), varargin{:});
    end

    function varargout = ehparaxial(beam, varargin)
      % Calculate the electric and magnetic field (in SI units)
      %
      % Usage
      %   [E, H, vswfData] = beam.ehparaxial(xy, ...)
      %
      % Returns
      %   - 3xNxM numeric array of values.  Where M is the number of beams.
      %
      % Parameters
      %   - xy (2xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of meters.  Packaged [x;y].
      %
      % See :class:`Beam` for additional information.

      [varargout{1:nargout}] = beam.doubleOutputHelper(...
          @(odata, varargin) odata.ehparaxial(varargin{:}), varargin{:});
    end
  end

  methods (Hidden)
    function [EH, vswfData] = singleOutputHelper(beam, ...
          target, xyz, varargin)
      % Helper for single output field calculation functions

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData());
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      target = @(b, varargin) target(...
          b.rotate(beam.rotation).translateXyz(beam.position), ...
          varargin{:});

      [EH, vswfData] = target(beam.data(1), xyz, unmatched{:});
      EH = repmat(EH, [1, 1, numel(beam.data)]);

      for ii = 1:numel(beam.data)
        [EH(:, :, ii), vswfData] = target(beam.data(ii), xyz, ...
            unmatched{:}, 'data', vswfData);
      end

      if strcmpi(beam.arrayType, 'coherent')
        EH = sum(EH, 3);
      end
    end

    function [E, H, vswfData] = doubleOutputHelper(beam, ...
          target, xyz, varargin)
      % Helper for double output field calculation functions

      p = inputParser;
      p.addParameter('data', ott.utils.VswfData());
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);
      
      target = @(b, varargin) target(...
          b.rotate(beam.rotation).translateXyz(beam.position), ...
          varargin{:});

      [E, H, vswfData] = target(beam.data(1), xyz, unmatched{:});
      E = repmat(E, [1, 1, numel(beam.data)]);
      H = repmat(H, [1, 1, numel(beam.data)]);

      for ii = 2:numel(beam.data)
        [E(:, :, ii), H(:, :, ii), vswfData] = target(beam.data(ii), ...
            xyz, unmatched{:}, 'data', vswfData);
      end

      if strcmpi(beam.arrayType, 'coherent')
        E = sum(E, 3);
        H = sum(H, 3);
      end
    end

    function val = defaultVisRangeInternal(beam)
      val = beam.data.defaultVisRange();
    end
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

      assert(numel(unique([val.index_medium])) == 1, ...
          'all beams must have the same refractive index');
      assert(numel(unique([val.omega])) == 1, ...
          'all beams must have the same optical frequency');

      % Store data
      beam.dataInternal = val;
    end
    function val = get.data(beam)
      val = beam.dataInternal;
    end

    function val = get.omega(beam)
      val = beam.data(1).omega;
    end

    function val = get.index_medium(beam)
      val = beam.data(1).index_medium;
    end
  end
end
