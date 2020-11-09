classdef Empty < ott.beam.Beam
% A beam with no fields.
% Inherits from :class:`Beam`.
%
% This class represents empty space with no fields.  This is useful
% for default parameters to functions or for unassigned elements of
% beam arrays (part of Hetrogeneous interface).
%
% Properties
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

    function [moment, ints, data] = intensityMoment(~, varargin)
      % Calculate moment of the beam intensity in the far-field.
      %
      % For :class:`Empty` beams, this returns zeros.
      %
      % Usage
      %   [moment, int, data] = beam.intensityMoment(...)

      data = [];
      moment = zeros(3, 1);
      ints = 0;
    end

    %
    % Field calculation methods (return zeros)
    %

    function [E, vswfData] = efield(~, xyz, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, vswfData] = beam.efield(xyz, ...)

      E = ott.utils.FieldVectorCart(zeros(size(xyz)));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfield(~, xyz, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfield(xyz, ...)

      H = ott.utils.FieldVectorCart(zeros(size(xyz)));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfield(~, xyz, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(xyz, ...)

      E = ott.utils.FieldVectorCart(zeros(size(xyz)));
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, vswfData] = efieldRtp(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, vswfData] = beam.efield(rtp, ...)

      E = ott.utils.FieldVectorCart(zeros(size(rtp)));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfieldRtp(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfield(rtp, ...)

      H = ott.utils.FieldVectorCart(zeros(size(rtp)));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfieldRtp(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(rtp, ...)

      E = ott.utils.FieldVectorCart(zeros(size(rtp)));
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, vswfData] = efarfield(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, vswfData] = beam.efarfield(xyz, ...)

      E = ott.utils.FieldVectorCart(zeros(3, size(rtp, 2)));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfarfield(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfarfield(xyz, ...)

      H = ott.utils.FieldVectorCart(zeros(3, size(rtp, 2)));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfarfield(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfarfield(xyz, ...)

      E = ott.utils.FieldVectorCart(zeros(3, size(rtp, 2)));
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, vswfData] = eparaxial(~, xy, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, vswfData] = beam.eparaxial(xyz, ...)

      E = ott.utils.FieldVectorCart(zeros(3, size(xy, 2)));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hparaxial(~, xy, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hparaxial(xyz, ...)

      H = ott.utils.FieldVectorCart(zeros(3, size(xy, 2)));
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehparaxial(~, xy, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(xyz, ...)

      E = ott.utils.FieldVectorCart(zeros(3, size(xy, 2)));
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end
  end
end

