classdef Scattered < ott.beam.Beam & ott.beam.abstract.Scattered
% Overload of the Beam class for scattered beams.
% Inherits from :class:`Beam` and :class:`abstract.Scattered`.
%
% This class provides simpler force and torque calculation methods
% which operate directly on the internal incident beam and ensure
% that the beam type is set appropriately before the calculation.
%
% For details on methods, see :class:`Beam`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods
    function varargout = force(ibeam, varargin)
      % Calculate change in linear momentum between beams.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   force = ibeam.force(other, ...)
      %
      % If other is not supplied, uses the incident beam for the calculation.

      p = inputParser;
      p.addOptional('other', [], @(x) isa(x, 'ott.beam.abstract.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if isempty(p.Results.other)
        if isempty(ibeam.incident_beam)
          error('Must supply other or have valid incident_beam');
        end

        [varargout{1:nargout}] = ibeam.incident_beam.force(...
            ibeam.total_beam, unmatched{:});
      else
        [varargout{1:nargout}] = force@ott.beam.Beam(ibeam, ...
            p.Results.other, unmatched{:});
      end
    end

    function varargout = torque(ibeam, varargin)
      % Calculate change in angular momentum between beams.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   torque = ibeam.torque(other, ...)
      %
      % If other is not supplied, uses the incident beam for the calculation.

      p = inputParser;
      p.addOptional('other', [], @(x) isa(x, 'ott.beam.abstract.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if isempty(p.Results.other)
        if isempty(ibeam.incident_beam)
          error('Must supply other or have valid incident_beam');
        end

        [varargout{1:nargout}] = ibeam.incident_beam.torque(...
            ibeam.total_beam, unmatched{:});
      else
        [varargout{1:nargout}] = torque@ott.beam.Beam(ibeam, ...
            p.Results.other, unmatched{:});
      end
    end

    function varargout = forcetorque(ibeam, varargin)
      % Calculate change in momentum between beams.
      %
      % Usage
      %   [force, torque] = ibeam.forcetorque(other, ...)
      %   Returns 3xN matrices for the force and torque in Cartesian
      %   coordinates.
      %
      % If other is not supplied, uses the incident beam for the calculation.
      %
      % Parameters
      %   - other (Beam|scat.Scatter) -- A beam to compare the force
      %     with or a particle with a scatter method.
      %
      %   - position (3xN numeric) -- Distance to translate beam before
      %     calculating the scattered beam using the T-matrix.
      %     Default: ``[]``.
      %
      %   - rotation (3x3N numeric) -- Angle to rotate beam before
      %     calculating the scattered beam using the T-matrix.
      %     Inverse rotation is applied to scattered beam, effectively
      %     rotating the particle.
      %     Default: ``[]``.

      p = inputParser;
      p.addOptional('other', [], @(x) isa(x, 'ott.beam.abstract.Beam'));
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      if isempty(p.Results.other)
        if isempty(ibeam.incident_beam)
          error('Must supply other or have valid incident_beam');
        end

        [varargout{1:nargout}] = ibeam.incident_beam.forcetorque(...
            ibeam.total_beam, unmatched{:});
      else
        [varargout{1:nargout}] = forcetorque@ott.beam.Beam(ibeam, ...
            p.Results.other, unmatched{:});
      end
    end
  end
end
