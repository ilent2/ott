classdef Shape < ott.scat.utils.ZeroScattered ...
    & ott.scat.utils.ShapeProperty ...
    & ott.scat.utils.HomogeneousRelative
% Shape-induced force approximation for homogeneous shapes.
% Inherits from :class:`ott.scat.utils.ZeroScattered`,
% :class:`ott.scat.utils.ShapeProperty` and
% :class:`ott.scat.utils.HomogeneousRelative`.
%
% The approximation involves calculating the electric field along the
% surface of the particle and estimating the force according to::
%
%     \langle f \rangle = \frac{1}{4} |E|^2 \Delta\epsilon \hat{n}
%
% where :math:`\hat{n}` is the normal to the surface,
% :math:`|E|^2` is the field amplitude and :math:`\Delta\epsilon` is
% the difference in permittivity of the particle and medium.
%
% To be accurate, the above should use the total field (incident+scatted),
% however, a zero-th order approximation is to use just the incident
% field, this approximation was shown to work rather well in
%
%   D. B. Phillips, et al., Shape-induced force fields in optical trapping
%   Nature Photonics volume 8, pages400–405(2014)
%   https://doi.org/10.1038/nphoton.2014.74
%
% Properties
%   - index_relative  -- Relative refractive index of particle.
%   - shape           -- Geometry associated with particle.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  methods
    function particle = Shape(shape, index_relative, varargin)
      % Construct a new homogeneous shape-induced force shape.
      %
      % Usage
      %   shape = Shape(shape, index_relative, ...)
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
  end

  methods (Hidden)
    function force = forceInternal(particle, beam, varargin)
      % Calculate the force using the shape-induced force approximation
      %
      % Implements::
      %
      %     \langle f \rangle = \frac{1}{4} |E|^2 \Delta\epsilon \hat{n}

      % Method can use any type with an efield method
      if ~isa(beam, 'ott.beam.Beam')
        beam = ott.beam.Beam(beam);
      end

      % Get surface locations
      [xyz, nxyz, dA] = particle.shape.surfPoints();

      % Calculate field at locations
      E = beam.efield(xyz);

      % Apply function for force calc to each beam in array
      epsilon = particle.index_relative.^2;
      func = @(x) sum(0.25 .* sum(abs(x).^2, 1) .* epsilon .* nxyz .* dA, 2);
      force = beam.arrayApply(func, E);

      % Combine incoherent force
      if isa(beam, 'ott.beam.utils.ArrayType')
        force = beam.combineIncoherentArray(force, 2);
        if iscell(force)
          force = cell2mat(force);
        end
      end
    end

    function t = torqueInternal(particle, beam, varargin)
      % Calculate torque

      % TODO: Implement this (same as geometric rays?)
      warning('Not yet implemented');
      t = [0;0;0];
    end

    function [f, t] = forcetorqueInternal(particle, beam, varargin)
      % Calculate force and torque

      f = particle.forceInternal(beam, varargin{:});
      t = particle.torqueInternal(beam, varargin{:});
    end
  end
end

