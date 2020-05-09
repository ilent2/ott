classdef ScatteredRay < ott.beam.Ray & ott.beam.Scattered
% Represents scattered ray objects
% Inherits from :class:`Ray` and :class:`abstract.Scattered`.
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
%
% Static methods
%   - empty             -- Construct an empty ScatteredRay array

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function beam = empty(type, varargin)
      % Construct an empty Scattered ray array
      %
      % Usage
      %   beam = ScatteredRay.empty(type, ...)
      %
      % Calls the class constructor with empty arrays.
      % Additional arguments are passed to the constructor.
      
      empt = zeros(3, 0);
      
      beam = ott.beam.ScatteredRay(type, ...
        'origin', empt, 'direction', empt, ...
        'polarisation', empt, 'field', empt(1, :), varargin{:});
    end
  end

  methods
    function beam = ScatteredRay(type, varargin)
      % Construct a new scattered ray representation
      %
      % Usage
      %   ray = ScatteredRay(type, ...)
      %
      % Parameters
      %   - incident_beam (Ray|empty) -- The incident beam or empty.
      %     Some functions require an incident beam.
      %     Default: ``[]``.
      %
      %   - type (enum) -- Type of scattered ray.
      %
      %   - like (Beam) -- Another Beam object to use for default
      %     properties.  If the object is a Ray object, copies the
      %     Ray properties (direction/field/...) too.
      %
      % All other parameters are passed to :class:`Ray`.
      
      p = inputParser;
      p.addParameter('incident_beam', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      beam = beam@ott.beam.Scattered(type, ...
          'incident_beam', p.Results.incident_beam);
      beam = beam@ott.beam.Ray(unmatched{:});
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
      
      % Setup hold
      washoldon = ishold();
      if ~washoldon
        hold('on');
      end

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
      
      % Return hold to previous state
      if ~washoldon
        hold('off');
      end

      if nargout > 0
        varargout{1} = h;
      end
    end
  end
end
