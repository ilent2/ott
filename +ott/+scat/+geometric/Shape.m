classdef Shape < ott.scat.utils.ShapeProperty ...
    & ott.scat.utils.RelativeMediumProperty
    & ott.scat.utils.BeamForce
% Describes scattering by a Shape instance.
% Inherits from :class:`ott.shapes.Shape`.
%
% This class describes scattering by Homogeneous shapes with a single
% refractive index value.
%
% Properties
%   - relativeMedium      -- Relative medium describing optical properties
%   - shape               -- Shape associated with particle
%
% Methods
%   - scatter             -- Calculate scattered rays

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods
    function particle = Shape(shape, relativeMedium, varargin)
      % Construct a new homogeneous geometric scattering shape.
      %
      % Usage
      %   particle = Shape(shape, relativeMedium, ...)
      %
      % Parameters
      %   - shape (ott.shapes.Shape) -- The geometric shape describing
      %     the shapes surface/boundary.
      %
      %   - index_relative (numeric) -- Relative refractive index of shape.

      p = inputParser;
      p.addOptional('shape', [], @(x) isa(x, 'ott.shapes.Shape'));
      p.addOptional('relativeMedium', [], ...
          @(x) isa(x, 'ott.beam.medium.Relative'));
      p.parse(varargin{:});

      particle.relativeMedium = p.Results.relativeMedium;
      particle.shape = p.Results.shape;
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

      inside = part.shape.insideXyz(ray.origin);
      assert(all(inside(:)) || all(~inside(:)), ...
        'Rays must be all inside or all outside shape');
      
      int_ray = ray;
      out_ray = ott.beam.ScatteredRay.empty('total', 'incident_beam', ray);
      
      empty = ott.beam.Ray.empty();
      if all(inside(:))
        int_ray = ott.beam.abstract.Array('array', [1, 2], empty, int_ray);
      else
        int_ray = ott.beam.abstract.Array('array', [1, 2], int_ray, empty);
      end

      for ii = 1:Niter

        new_int_ray = ott.beam.abstract.Array('array', size(int_ray));
        for jj = 1:2
          
          % Calculate intersections
          [locs, norms] = part.shape.intersect(int_ray(jj));

          % Split rays into interacting and non-interacting
          no_intersect = any(isnan(locs), 1);
          if any(no_intersect)
            out_ray(end+(1:sum(no_intersect))) = int_ray.beams{jj}(no_intersect); %#ok<AGROW>
            int_ray.beams{jj} = int_ray.beams{jj}(~no_intersect);

            locs(:, no_intersect) = [];
            norms(:, no_intersect) = [];
          end
          
          if isempty(int_ray(jj))
            continue;
          end
          
          % Construct planes for interacting rays
          plane = ott.scat.geometric.Plane(norms, part.index_relative, ...
              'position', locs);

          % Calculate scattered ray directions
          [rplane, tplane] = plane.scatter(int_ray(jj));
          
          % Form new ray object
          if jj == 1
            % Rays were exterior, rplane -> exterior(1), tplane -> interior(2)
            new_int_ray(1) = rplane;
            new_int_ray(2) = tplane;
          else
            % Rays were interior, tplane -> exterior(1), rplane -> interior(2)
            new_int_ray(1) = [new_int_ray(1), tplane];
            new_int_ray(2) = [new_int_ray(2), rplane];
          end
        end
        
        int_ray = new_int_ray;
        
        % Prune rays bellow threshold power
        if ~isempty(p.Results.prune_threshold)
          for jj = 1:numel(int_ray)
            discard = int_ray(jj).intensity < p.Results.prune_threshold;
            int_ray.beams{jj}(discard) = [];
          end
        end

        % Check we have any rays left
        if isempty(new_int_ray(1)) && isempty(new_int_ray(2))
          break;
        end

        % Check early stopping threshold
        if sum([int_ray.power{:}]) < p.Results.stop_threshold
          int_ray.beams = {};
          break;
        end
      end
    end
  end
end
