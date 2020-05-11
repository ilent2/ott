classdef Plane < ott.shapes.Plane & ott.scat.utils.Particle ...
    & ott.scat.utils.HomogeneousRelative ...
    & ott.scat.utils.BeamForce
% Describes how a infinite plane scatters a plane wave.
% Inherits from :class:`ott.shapes.Plane`.
%
% The refractive index of the medium above the surface (in the position
% direction of the normal vector) is related to the index inside the
% surface (negative to the noraml) by:
%
%   index_relative = (negative) ./ (posative)
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
    function plane = Plane(normal, index_relative, varargin)
      % Construct a new plane instance
      %
      % Usage
      %   plane = Plane(normal, index_relative, ...)
      %
      % Optional named arguments are passed to :class:`ott.shapes.Plane`.
      
      plane = plane@ott.shapes.Plane(normal, varargin{:});
      plane = plane@ott.scat.utils.HomogeneousRelative(index_relative);
    end

    function [rbeam, tbeam] = scatter(plane, beam)
      % Calculate reflected and transmitted beams
      
      import ott.utils.cross;
      import ott.utils.dot;
      
      % Handle arrays of beams
      % TODO: We may be able to do this more efficiently by creating
      % a array with pre-scaled values
      if isa(beam, 'ott.beam.abstract.Array')
        rbeam = ott.beam.Array(beam.array_type, size(beam));
        tbeam = ott.beam.Array(beam.array_type, size(beam));
        for ii = 1:numel(beam)
          [rbeam(ii), tbeam(ii)] = plane.scatter(beam(ii));
        end
        return;
      end

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
      
      % Get the relative index (flipping if inside)
      index_relative = plane.index_relative;
      mask = sign(dot(plane.normal, ki)) > 0;
      if isscalar(index_relative)
        index_relative = repmat(index_relative, size(mask));
      end
      index_relative(mask) = 1.0./index_relative(mask);

      % Calculate transmitted wave-vector (direction)
      % Orthogonal components are unchanged
      % Normal component is scaled by relative index
      ko = ki - dot(plane.normal, ki).*plane.normal;  % orthogonal
      kn2 = index_relative.^2 .* dot(ki, ki) - dot(ko, ko); % normal
      kt = ko + sqrt(kn2) .* plane.normal .* sign(dot(plane.normal, ki));

      % Calculate polarisation
      polarisation = double(beam.polarisation);
      
      % Calculate vector perpendicular to plane
      % If ki normal to plane, use polarisation direction
      svec = cross(beam.direction, plane.normal);
      tol = 1.0e-6;
      svecmask = vecnorm(svec) < tol;
      svec(:, svecmask) = polarisation(:, svecmask);
      svec = svec ./ vecnorm(svec);

      % Split the field (polarisation) into s and p vectors
      Es = dot(polarisation, svec);
      Ep = vecnorm(polarisation - Es .* svec);
      
      % Multiple by the field value
      if size(beam.field, 1) == 1
        Es = Es .* beam.field;
        Ep = Ep .* beam.field;
      else
        sz = size(beam.field);
        Es = Es .* reshape(beam.field(1, :), [1, sz(2:end)]);
        Ep = Ep .* reshape(beam.field(1, :), [1, sz(2:end)]);
        
        % Do the same again for the orthogonal polarisation
        orth_pol = cross(polarisation, beam.direction);
        Es2 = dot(orth_pol, svec);
        Ep2 = vecnorm(orth_pol - Es .* svec);
        Es = Es + Es2 .* reshape(beam.field(2, :), [1, sz(2:end)]);
        Ep = Ep + Ep2 .* reshape(beam.field(2, :), [1, sz(2:end)]);
      end

      % Calculate refractive index
      n1 = beam.medium.index;
      n2 = beam.medium.index .* index_relative;

      % Calculate the Fresnel coefficients
      kix = dot(plane.normal, ki);
      ktx = dot(plane.normal, kt);
      [Sr, St] = plane.fresnelS(kix, ktx, n1, n2);
      [Pr, Pt] = plane.fresnelP(kix, ktx, n1, n2);

      % Calculate transmitted and reflected pvec directions
      pvecr = cross(svec, kr)./vecnorm(kr);
      pvect = cross(svec, kt)./vecnorm(kt);

      % Generate the reflected and transmitted vectors
      % TODO: Not sure if real(kt) is the right thing to do???
      rbeam = ott.beam.abstract.PlaneWave('direction', kr./vecnorm(kr), ...
          'polarisation', pvecr, ...
          'field', [Pr .* Ep; Sr .* Es], ...
          'origin', beam.origin, ...
          'like', beam);
      tbeam = ott.beam.abstract.PlaneWave('direction', real(kt)./vecnorm(real(kt)), ...
          'polarisation', pvect, ...
          'field', [Pt .* Ep; St .* Es], ...
          'index', unique(n2), ...
          'origin', beam.origin, ...
          'like', beam);
    end
  end
end
