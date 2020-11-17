classdef Empty < ott.beam.Beam & ott.beam.properties.IndexOmegaProps
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
    function bm = Empty(varargin)
      % Construct a new empty beam instance.
      %
      % Usage
      %   beam = ott.beam.Empty(...)
      %
      % All parameters passed to :class:`Beam` constructor.

      [omega, index, unmatched] = ott.beam.properties. ...
          IndexOmegaProps.parseArgs(varargin{:});

      bm = bm@ott.beam.properties.IndexOmegaProps(omega, index);
      bm = bm@ott.beam.Beam(unmatched{:});
    end
    
    function bsc = ott.bsc.Bsc(~, varargin)
      % Construct empty Bsc instance
      bsc = ott.bsc.Bsc();
    end

    %
    % Force calculation methods
    %
    
    function varargout = force(ebeam, other, varargin)
      % Calculate force, defering to other beam force method
      %
      % Usage
      %   [F, ...] = emptyBeam.force(other, ...)
      
      [varargout{1:nargout}] = other.force(ebeam, varargin{:});
      varargout{1} = -varargout{1};
    end
    
    function varargout = torque(ebeam, other, varargin)
      % Calculate torque, defering to other beam torque method
      %
      % Usage
      %   [F, ...] = emptyBeam.torque(other, ...)
      
      [varargout{1:nargout}] = other.torque(ebeam, varargin{:});
      varargout{1} = -varargout{1};
    end
    
    function varargout = spin(ebeam, other, varargin)
      % Calculate spin, defering to other beam spin method
      %
      % Usage
      %   [F, ...] = emptyBeam.spin(other, ...)
      
      [varargout{1:nargout}] = other.spin(ebeam, varargin{:});
      varargout{1} = -varargout{1};
    end

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

      E = ott.utils.FieldVectorCart(zeros(size(xyz)), xyz);
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfield(~, xyz, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfield(xyz, ...)

      H = ott.utils.FieldVectorCart(zeros(size(xyz)), xyz);
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfield(~, xyz, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(xyz, ...)

      E = ott.utils.FieldVectorCart(zeros(size(xyz)), xyz);
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

      E = ott.utils.FieldVectorSph(zeros(size(rtp)), rtp);
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfieldRtp(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfield(rtp, ...)

      H = ott.utils.FieldVectorSph(zeros(size(rtp)), rtp);
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfieldRtp(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfield(rtp, ...)

      E = ott.utils.FieldVectorSph(zeros(size(rtp)), rtp);
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

      % Ensure rtp size is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      E = ott.utils.FieldVectorSph(zeros(3, size(rtp, 2)), rtp);
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [H, vswfData] = hfarfield(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [H, vswfData] = beam.hfarfield(xyz, ...)

      % Ensure rtp size is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      H = ott.utils.FieldVectorSph(zeros(3, size(rtp, 2)), rtp);
      if nargout == 2
        vswfData = ott.utils.VswfData();
      end
    end

    function [E, H, vswfData] = ehfarfield(~, rtp, varargin)
      % Returns zeros the same size as input coordinates.
      %
      % Usage
      %   [E, H, vswfData] = beam.ehfarfield(xyz, ...)

      % Ensure rtp size is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      E = ott.utils.FieldVectorSph(zeros(3, size(rtp, 2)), rtp);
      H = E;
      if nargout == 3
        vswfData = ott.utils.VswfData();
      end
    end
  end
end

