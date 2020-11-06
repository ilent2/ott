classdef Empty < ott.beam.Beam
% A beam with no fields.
% Inherits from :class:`Beam`.
%
% This class represents empty space with no fields.  This is useful
% for default parameters to functions or for unassigned elements of
% beam arrays (part of Hetrogeneous interface).
%
% Properties
%   - data      - An empty beam array
%   - position  - (3x1 numeric) Position of the empty space.
%   - rotation  - (3x3 numeric) Orientation of the empty space.
%   - index_medium    -- Refractive index of the medium
%   - wavelength      -- Wavelength in medium [m]
%   - wavenumber      -- Wavenumber in medium [1/m]
%   - omega           -- Optical angular frequency of light [1/s]
%   - speed           -- Speed of light in the medium [m/s]
%   - speed0          -- Speed of light in vacuum [m/s]
%
% Methods
%   - efield, hfield, ehfield, ...    -- Returns zeros

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    data = ott.beam.Beam.empty(0);
  end

  methods
    function beam = Empty(varargin)
      % Construct a new empty beam instance.
      %
      % Usage
      %   beam = ott.beam.Empty(...)
      %
      % All parameters passed to :class:`Beam` constructor.

      beam = beam@ott.beam.Beam(varargin{:});
    end

    %
    % Force calculation methods
    %

    function [moment, ints, data] = intensityMoment(beam, varargin)
      % Calculate moment of the beam intensity in the far-field.
      %
      % For :class:`Empty` beams, this returns zeros.
      %
      % Usage
      %   [moment, int, data] = beam.intensityMoment(...)

      data = [];
      moment = zeros(3, 1);
      int = 0;
    end

    %
    % Field calculation methods (return zeros)
    %

    function [E, vswfData] = efield(beam, xyz, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, vswfData] = beam.efield(xyz, ...)

      E = zeros(size(xyz));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfield(beam, xyz, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfield(xyz, ...)

      H = zeros(size(xyz));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfield(beam, xyz, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(xyz, ...)

      E = zeros(size(xyz));
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, vswfData] = efieldRtp(beam, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, vswfData] = beam.efield(rtp, ...)

      E = zeros(size(rtp));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfieldRtp(beam, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfield(rtp, ...)

      H = zeros(size(rtp));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfieldRtp(beam, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(rtp, ...)

      E = zeros(size(rtp));
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, vswfData] = efarfield(beam, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, vswfData] = beam.efarfield(xyz, ...)

      E = zeros(3, size(rtp, 2));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfarfield(beam, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfarfield(xyz, ...)

      H = zeros(3, size(rtp, 2));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfarfield(beam, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfarfield(xyz, ...)

      E = zeros(3, size(rtp, 2));
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, vswfData] = eparaxial(beam, xy, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, vswfData] = beam.eparaxial(xyz, ...)

      E = zeros(3, size(xy, 2));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hparaxial(beam, xy, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hparaxial(xyz, ...)

      H = zeros(3, size(xy, 2));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehparaxial(beam, xy, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(xyz, ...)

      E = zeros(3, size(xy, 2));
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end
  end

  methods % Getters/setters
    function beam = set.data(beam, ~)
      warning('Setting beam data has no effect for Empty beams');
    end
  end
end

