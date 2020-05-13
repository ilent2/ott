classdef Dipole < ott.scat.utils.Particle & ott.utils.RotationPositionProp
% Single point dipole approximation for the force on a particle.
% Inherits from :class:`ott.scat.utils.Particle`.
%
% When the particle is small compared to the wavelength, the fields
% across the particle can be assumed to be uniform, i.e. the particle
% is in the Rayleigh scattering regime and can be modelled as a small
% polarizable dipole.  In this case, the dipole polarisation is given by::
%
%   p = \alpha \vec{E}
%
% where :math:`\alpha` is the dipole polarizability and :math:`\vec{E}`
% is the electric field.
%
% The forces are given by::
%
%   F_{scat} = \frac{n}{c} C_{pr} \vec{I}
%
% and::
%
%   F_{grad} = \alpha \nabla\vec{I}
%
% and :math:`\vec{I}` is the time average of the Poynting vector,
% which for paraxial beams often simplifies to :math:`\hat{z} I`,
% :math:`n` is the medium refractive index, :math;`c` is the
% speed of light in vacuum, :math:`C_{pr}` is the radiation pressure
% cross-section of the particle.
%
% The dipole limit is effectively the :math:`N_{max} = 1` limit of
% the Mie scattering description.  For further details see
%
%   Nieminen et al., Approximate and exact modeling of optical trapping
%   Proc. SPIE (7762), Optical Trapping and Optical Micromanipulation VII,
%   77622V, (2010), https://doi.org/10.1117/12.861880
%
% and
%
%   Yasuhiro Harada and Toshimitsu Asakura,
%   Optics Communications, Volume 124, Issues 5–6, 1996, Pages 529-541
%   https://doi.org/10.1016/0030-4018(95)00753-9
%
% Methods
%   - force       -- Calculate the force acting on the dipole
%   - force_scat  -- Calculate the scattered force only
%   - force_grad  -- Calculate the gradient force only
%   - scatter     -- Generates a :class:`ScatteredBeam`
%
% Properties
%   - polarizability    -- Polarizability (scalar | 3x1 | 3x3 numeric)
%   - rpcrosssection    -- Radiation pressure cross-section (1|3x1|3x3)
%
% See also Dielectric

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    polarizability      % Polarizability (scalar | 3x1 | 3x3 numeric)
    rpcrosssection      % Radiation pressure cross-section (1|3x1|3x3)
  end

  methods (Static)
    function val = radiativeReaction(alp)

      if isvector(alp)
        alp = diag(alp);
      end

      index_relative = 1.2;
      warning('N-rel and k hard coaded for now, fit it!');

%       n = 1.0;
      epsU = index_relative.^2;
      a = 0.5;
      k = 2*pi/1.0;

      % Equation 2.03 from Draine 1988.
      % The Discrete-Dipole approximation and its application to
      % interstellar graphite grains.  The astrophysical journal,
      % 333 (848--872), 1988.
      %
      % TODO: I don't think this is the right equation
%       val = alp .* inv(1 - 2i./(4*pi*n) .* k.^3 .* (epsU - 1)./(epsU + 2));

      % This one is Equation 11, Yasuhiro Harada and Toshimitsu Asakura
      val = 8/3 .* pi .* (k .* a).^4 .* a.^2 .* ((epsU - 1)./(epsU + 2)).^2;

      warning('Might not have the right cpr');

    end
  end

  methods
    function dipole = Dipole(polarizability, varargin)
      % Construct a new polarizable dipole instance
      %
      % Usage
      %   dipole = Dipole(polarizability, ...)
      %
      % Parameters
      %   - polarizability (numeric) -- Dipole polarizability.
      %     Can either be a scalar, 3x1 vector or 3x3 matrix.
      %
      % Optional named arguments
      %   - rpcrosssection (numeric) -- Radiation pressure
      %     cross section.  Can be scalar, 3x1 or 3x3.
      %     Default: ``Dipole.radiativeReaction(polarizability)``.

      % Parse arguments
      p = inputParser;
      p.addParameter('rpcrosssection', ...
        ott.scat.dipole.Dipole.radiativeReaction(polarizability));
      p.parse(varargin{:});

      % Store properties
      dipole.polarizability = polarizability;
      dipole.rpcrosssection = p.Results.rpcrosssection;
    end
  end

  methods (Hidden)

    function sbeam = scatterInternal(dipole, beam, xyz)
      % Generates a :class:`ScatteredBeam`
      %
      % Usage
      %   sbeam = dipole.scatter(beam, xyz)
      %
      % Parameters
      %   - beam -- A :class:`ott.optics.beam.Beam`
      %   - xyz (3x1 numeric) -- Dipole position in beam

      % Defer all argument checks to ScatteredBeam
      % TODO: Should we have a radiating Dipole beam type?
      sbeam = ott.scat.dipole.ScatteredBeam(beam, dipole, xyz);
    end

    function f = forceInternal(dipole, beam)
      % Calculate the total force on a dipole in a beam

      % Cast beam to a real beam
      if ~isa(beam, 'ott.beam.Beam')
        beam = ott.beam.Beam(beam);
      end

      f = dipole.force_scat(beam) + dipole.force_grad(beam);
    end

    function t = torqueInternal(dipole, beam)

      t = zeros(3, 1);
      warning('Not yet implemented');

    end

    function f = force_scat(dipole, beam)
      % Calculate the scattered force on a dipole in a beam
      %
      % Usage
      %   f = dipole.force_scat(beam)
      %
      % Parameters
      %   - beam -- A :class:`ott.optics.beam.Beam`

      % Calculate Poynting vector
      S = beam.poynting(dipole.position);

      % Calculate force (Harada Eq 10)
      f = 1./beam.medium.speed .* dipole.rpcrosssection .* 0.5.*real(S);
    end

    function f = force_grad(dipole, beam)
      % Calculate the gradient force on a dipole in a beam
      %
      % Usage
      %   f = dipole.force_grad(beam)
      %
      % Parameters
      %   - beam -- A :class:`ott.optics.beam.Beam`

      % Calculate electric field at positions
      E = beam.efield(dipole.position);

      % Calculate the E-field gradient
      gradE = beam.egradient(dipole.position);

      % Calculate polarisation (Harada Eq 9)
      % TODO: Should this be the total or partial polarizability
      if isvector(dipole.polarizability) || isscalar(dipole.polarizability)
        p = dipole.polarizability .* E;
      else
        % Polarizability is 3x3 matrix
        p = dipole.polarizability * E;
      end

      % Calculate force (Harada Eq 13)
      p = reshape(p, 3, 1, []);
      f = 0.5 .* real(reshape(sum(p .* gradE, 1), 3, []));

    end
  end
end

