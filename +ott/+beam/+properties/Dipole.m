classdef Dipole < ott.beam.properties.Material
% Properties of dipoles.
% Inherits from :class:`ott.beam.properties.Beam`.
%
% Properties
%   - position      -- Position of the dipole collection
%   - rotation      -- Orientation of the dipole collection
%
% Methods
%   - setDipoles    -- Set dipole positions/polarizations
%
% Abstract properties
%   - location      -- Location(s) of individual dipole(s)
%   - polarization  -- Polarization of the dipole(s)

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Abstract, SetAccess=protected)
    location      % Location(s) of individual dipole(s)
    polarization  % Polarization of the dipole(s)
  end

  methods (Static)
    function args = likeProperties(other, args)
      % Add like-properties to argument list
      if isa(other, 'ott.beam.properties.Dipole')
        args = ott.utils.addDefaultParameter(...
            'location', other.location, args);
        args = ott.utils.addDefaultParameter(...
            'polarization', other.polarization, args);
      end
      args = ott.beam.properties.Material.likeProperties(other, args);
    end
  end

  methods
    function beam = Dipole(varargin)
      % Construct dipole parameters
      %
      % Usage
      %   beam = Dipole(location, polarization, ...)
      %   Parameters can also be passed as named arguments

      p = inputParser;
      p.addOptional('location', [], @isnumeric);
      p.addOptional('polarization', [], @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.properties.Material(unmatched{:});
      beam = beam.setDipoles(p.Results.location, p.Results.polarization);
    end

    function beam = setDipoles(beam, location, polarization)
      % Set the dipole position and polarization data
      %
      % Usage
      %   beam = beam.setDipoles(location, polarization)
      %
      % Parameters
      %   - location (3xN numeric) -- Locations of dipoles
      %   - polarization (3NxM) -- Dipole polarizations sorted
      %     packaged in [x1;y1;z1; x2;y2;z2; ...] order.
      %
      % Check implementation for supported array sizes.

      ott.utils.nargoutCheck(beam, nargout);

      assert(isnumeric(location) && ismatrix(location) && size(location, 1) == 3, ...
        'location must be 3xN numeric matrix');
      assert(isnumeric(polarization) && ismatrix(polarization) ...
          && size(polarization, 1) == numel(location), ...
          'polarization must be 3NxM numeric matrix');

      beam.location = location;
      beam.polarization = polarization;
    end
  end
end
