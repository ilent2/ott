classdef (Abstract) Beam
% Base class for beam approximations
%
% Methods
%   - efield    -- Calculate electric field
%   - hfield    -- Calculate magnetic field
%   - ehfield   -- Calculate electric and magnetic fields
%   - visualise -- Generate a visualisation of the fields
%
% Properties
%   - power     -- The power of the beam (may be infinite)
%
% Abstract methods
%   - efieldInternal    -- Called by efield
%   - hfieldInternal    -- Called by hfield
%   - getBeamPower      -- get method called by dependent property power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties (Dependent)
    power           % The power of the beam (may be infinite)
  end

  methods (Abstract)
    efieldInternal      % Method called by efield
    hfieldInternal      % Method called by hfield
    getBeamPower        % Get the beam power
  end

  methods
    function E = efield(beam, xyz)
      % Calculate E and H field
      %
      % Usage
      %   E = beam.efield(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).

      E = beam.hfieldInternal(xyz);
    end

    function H = hfield(beam, xyz)
      % Calculate E and H field
      %
      % Usage
      %   H = beam.hfield(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).

      H = beam.efieldInternal(xyz);
    end

    function [E, H] = ehfield(beam, xyz)
      % Calculate E and H field
      %
      % Usage
      %   [E, H] = beam.ehfield(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).

      % Default implementation just calls each part, but for some
      % beams it may be more efficient to calculate these together,
      % in which case this method should be over-written.
      E = beam.efield(xyz);
      H = beam.hfield(xyz);
    end

    function visualise(beam, varargin)
      % Generate a visualisation of the beam
      %
      % TODO: This will be fleshed out in future releases
      % including moving the field data function from Bsc into
      % this class and making it a static method.
      %
      % Maybe also move other parts of ott.Bsc.visualise here too?

      % Generate grid of points
      x = linspace(-1, 1, 100).*beam.wavelength*10;
      y = linspace(-1, 1, 100).*beam.wavelength*10;
      [xx, yy, zz] = meshgrid(x, 0, y);
      %[xx, yy, zz] = meshgrid(x, y, 0);
      xyz = [xx(:), yy(:), zz(:)].';

      % Calculate field
      E = sum(abs(beam.efield(xyz)).^2, 1);

      % Draw image
      imagesc(x, y, squeeze(reshape(E(1, :), size(xx))));
    end
  end

  methods (Hidden)
    function beam = setBeamPower(beam, val)
      % Function to set the beam power (if supported)
      % Override this function if your beam supports this feature
      error('Setting beam power not supported');
    end
  end

  methods
    function val = get.power(beam)
      val = beam.getBeamPower();
    end
    function beam = set.power(beam, val)
      beam = beam.setBeamPower(val);
    end
  end
end

