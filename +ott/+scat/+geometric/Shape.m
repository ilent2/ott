classdef Shape < ott.scat.utils.ShapeProperty ...
    & ott.scat.utils.HomogeneousRelative ...
    & ott.scat.utils.BeamForce
% Describes scattering by a Shape instance.
% Inherits from :class:`ott.shapes.Shape`.
%
% This class describes scattering by Homogeneous shapes with a single
% refractive index value.
%
% Properties
%   - index_relative      -- Refractive index of particle relative to medium
%   - shape               -- Shape associated with particle
%
% Methods
%   - scatter             -- Calculate scattered rays

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods
    function particle = Shape(shape, index_relative, varargin)
      % Construct a new homogeneous geometric scattering shape.
      %
      % Usage
      %   particle = Shape(shape, index_relative, ...)
      %
      % Parameters
      %   - shape (ott.shapes.Shape) -- The geometric shape describing
      %     the shapes surface/boundary.
      %
      %   - index_relative (numeric) -- Relative refractive index of shape.

      p = inputParser;
      p.parse(varargin{:});

      particle = particle@ott.scat.utils.HomogeneousRelative(index_relative);
      particle = particle@ott.scat.utils.ShapeProperty(shape);
    end

    function [out_ray, int_ray] = scatter(part, ray, varargin)
      % Calculate scattering from an incident ray
      %
      % Usage
      %   [out_ray, int_ray] = particle.scatter(ray)
      %
      % Parameters
      %   - ray (ott.beam.Ray) -- Incident ray set
      %
      % Outputs
      %   - out_ray (ott.beam.ScatteredRay) -- Outgoing scattered rays
      %     and rays not scattered by the particle.
      %     Rays leaving at the same time are placed in the same row.
      %
      %   - int_ray (ott.beam.ScatteredRay) -- Rays still interacting
      %     with the particle after iteration completes.
      %
      % Optional named arguments
      %   - max_iterations (numeric) -- Maximum number of iterations.
      %     Default: ``10``.
      %
      %   - stop_threshold (numeric) -- Threshold to stop iterations.
      %     Default: ``1e-6``.
      %
      %   - prune_threshold (numeric) -- Threshold to use when discarding
      %     rays from further scattering calculations.
      %     Default: ``1e-6``.

      p = inputParser;
      p.addParameter('max_iterations', 10);
      p.addParameter('stop_threshold', 1e-6);
      p.addParameter('prune_threshold', 1e-6);
      p.parse(varargin{:});

      Niter = p.Results.max_iterations;
      assert(Niter > 0 && isscalar(Niter) && round(Niter) == Niter, ...
          'max_iterations must be positive integer');
        
      if ~isa(ray, 'ott.beam.Ray')
        ray = ott.beam.Ray(ray);
      end

      int_ray = ray;
      out_ray = ott.beam.ScatteredRay.empty(ray, 'total');

      for ii = 1:Niter

        % Calculate intersections
        [locs, norms] = part.shape.intersect(int_ray);

        % Split rays into interacting and non-interacting
        no_intersect = any(isnan(locs), 1);
        if any(no_intersect)
          out_ray(end+(1:sum(no_intersect))) = int_ray(no_intersect); %#ok<AGROW>
          int_ray = int_ray(~no_intersect);
          
          locs(:, no_intersect) = [];
          norms(:, no_intersect) = [];
        end

        % Check we have work to do
        if isempty(int_ray)
          break;
        end

        % Construct planes for interacting rays
        plane = ott.scat.planewave.Plane(norms, part.index_relative, ...
            'position', locs);

        % Calculate scattered ray directions
        [rplane, tplane] = plane.scatter(ray);

        % Form new ray object
        int_ray = ott.beam.Ray([rplane, tplane]);

        % Prune rays bellow threshold power
        if ~isempty(p.Results.prune_threshold)
          discard = int_ray.intensity < p.Results.prune_threshold;
          int_ray(discard) = [];
        end

        % Check early stopping threshold
        if int_ray.power < p.Results.stop_threshold
          int_ray(:) = [];
          break;
        end
      end
    end
  end
end
