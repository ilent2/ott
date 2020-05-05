classdef Plane < ott.shapes.Plane
% Describes how a infinite plane scatters a plane wave.
% Inherits from :class:`ott.shapes.Plane`.
%
% Properties
%   - normal          -- Vector normal to the plane surface
%   - index_relative  -- Relative refractive index of plane to medium
%
% Methods
%   - scatter         -- Calculate scattered plane wave beams
%
% Static methods
%   - fresnelS        -- Calculate S-direction Fresnel coefficient
%   - fresnelP        -- Calculate P-direction Fresnel coefficient

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file


% TODO: Ensure testScatterIndexMatched runs
% TODO: Inherit from another class for index_relative and geometry
% TODO: Move somewhere else (this isn't geometric optics) (where?)
% TODO: Inherit a geometric optics class from this class (how? why?)
% TODO: Other particles?  Other shapes?
% TODO: Implement a abstract plane wave class (use that instead? why?)
%   * Do we need an analytic beam and analytic scat sub-package?
%   * Mie is a analytic method for spheres, so why put it in T-matrix?
%   * Maybe we should have a ``+planewave`` sub-package for plane wave
%     basis scattering and solutions (different from geometric)
%   * Still, what about beams?  We also have annuarl, Top-hat and Davis?

  properties
    index_relative  % Relative refractive index
  end

  methods (Static)
    function [Sr, St] = fresnelS(kix, ktx, n1, n2)
      % Equations 33.51 and 33.52 from Feynman
      Sr = (kix - ktx)./(kix + ktx);
      St = 2.*kix ./ (kix + ktx);
    end

    function [Pr, Pt] = fresnelP(kix, ktx, n1, n2)
      % Equations 33.53 and 33.54 from Feynman
      Pr = (n2.^2 .* kix - n1.^2 .* ktx)./(n2.^2 .* kix + n1.^2 .* ktx);
      Pt = 2.*n1.*n2.*kix ./ (n2.^2 .* kix + n1.^2 .* ktx);
    end
  end

  methods
    function plane = Plane(index_relative, normal)
      % Construct a new plane instance
      plane.index_relative = index_relative;
      plane.normal = normal;
    end

    function [rbeam, tbeam] = scatter(plane, beam)
      % Calculate reflected and transmitted beams

      % Cast the beam to a plane wave
      if ~isa(beam, 'ott.beam.abstract.PlaneWave')
        beam = ott.beam.abstract.PlaneWave(beam);
      end

      % Get incident wave-vector (direction)
      ki = beam.wavevector;

      % Calculate reflected wave-vector (direction)
      % wavenumber is the same since the medium is the same
      % Direction of the normal component changes (reflected)
      kr = ki - 2.*dot(plane.normal, ki).*plane.normal;

      % Calculate transmitted wave-vector (direction)
      % Orthogonal components are unchanged
      % Normal component is scaled by relative index
      ko = ki - dot(plane.normal, ki).*plane.normal;  % orthogonal
      kn2 = plane.index_relative.^2 .* dot(ki, ki) - dot(ko, ko); % normal
      kt = ko + sqrt(kn2) .* plane.normal;

      % Calculate polarisation (phase-shifted by distance from beam
      % origin to plane)
      polarisation = beam.polarisation;

      % TODO: Use field vector

      % Split the field (polarisation) into s and p vectors
      svec = cross(beam.direction, plane.normal);
      Es = dot(polarisation, svec);
      Ep = vecnorm(polarisation - Es .* svec);

      n1 = beam.index_medium;
      n2 = beam.index_medium .* plane.index_relative;

      % Calculate the Fresnel coefficients
      kix = dot(plane.normal, ki);
      ktx = dot(plane.normal, kt);
      [Sr, St] = plane.fresnelS(kix, ktx, n1, n2);
      [Pr, Pt] = plane.fresnelP(kix, ktx, n1, n2);

      % Calculate transmitted and reflected pvec directions
      pvecr = cross(svec, kr)./vecnorm(kr);
      pvect = cross(svec, kt)./vecnorm(kt);

      % Generate the reflected and transmitted vectors
      rbeam = ott.beam.abstract.PlaneWave('direction', kr./vecnorm(kr), ...
          'polarisation', Sr .* Es .* svec + Pr .* Ep .* pvecr, ...
          'index_medium', n1, ...
          'origin', beam.origin);
      tbeam = ott.beam.abstract.PlaneWave('direction', kt./vecnorm(kt), ...
          'polarisation', St .* Es .* svec + Pt .* Ep .* pvect, ...
          'index_medium', n2, ...
          'origin', beam.origin);
    end
  end
end
