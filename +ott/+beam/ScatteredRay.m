classdef ScatteredRay < ott.beam.Ray
% Represents scattered ray objects
% Inherits from :class:`Ray`.
%
% This class contains information about the incident and scattered ray.
% This makes it easier to produce visualisations.
%
% Properties
%   - incident_beam     -- The incident ray object (can be set to [])
%
% Methods
%   - visualiseRays     -- Plot only scattered beam rays
%   - visualiseAllRays  -- Plot incident beam and scattered beam rays

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    incident_beam       % The incident ray object (can be set to [])
  end

  methods
    function beam = ScatteredRay(incident_beam, varargin)
      % Construct a new scattered ray representation
      %
      % Usage
      %   ray = ScatteredRay(incident_beam, ...)
      %
      % Parameters
      %   - incident_beam (Ray|empty) -- The incident beam or empty.
      %     Some functions require an incident beam.
      %
      % All other parameters are passed to :class:`Ray`.

      beam = beam@ott.beam.Ray(varargin{:});
      beam.incident_beam = incident_beam;
    end

    function varargout = visualiseAllRays(beam, varargin)
      % Generate a visualisation of incident beams and scattered rays
      %
      % If the incident beam is a Scattered ray, also calls the
      % visualiseAllRays on the incident beam.
      %
      % Usage
      %   h = beam.visualiseRays(...)
      %   Returns the quiver handles.
      %
      % Other arguments are passed to visualiseRays.

      % Visualise our rays
      h = beam.visualiseRays(varargin{:});

      % Loop over remaining beams
      last_beam = beam;
      next_beam = beam.incident_beam;
      while ~isempty(next_beam)

        % Defer to visualiseRays for visualisation
        ray_lengths = vecnorm(last_beam.origin - next_beam.origin);
        h = [h, next_beam.visualiseRays(varargin{:}, ...
            'ray_lengths', ray_lengths)];

        last_beam = next_beam;
        next_beam = [];
        if isa(last_beam, 'ott.beam.ScatteredRay')
          next_beam = last_beam.incident_beam;
        end
      end

      if nargout == 1
        varargout{1} = h;
      end
    end
  end

  methods % Getters/setters
    function beam = set.incident_beam(beam, val)
      assert(isempty(val) || isa(val, 'ott.beam.Ray'), ...
        'Incident beam must be empty or a Ray object');
      beam.incident_beam = val;
    end
  end
end
