classdef Dipole
% Single point dipole approximation for the force on a particle
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
      % TODO: Document where this is from

      if isvector(alp)
        alp = diag(alp);
      end

       % TODO: Implement
      error('Not yet implemented');

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
        ott.optics.dipole.Dipole.radiativeReaction(polarizability));
      p.parse(varargin{:});

      % Store properties
      dipole.polarizability = polarizability;
      dipole.rpcrosssection = p.Results.rpcrosssection;
    end

    function sbeam = scatter(dipole, beam, xyz)
      % Generates a :class:`ScatteredBeam`
      %
      % Usage
      %   sbeam = dipole.scatter(beam, xyz)
      %
      % Parameters
      %   - beam -- A :class:`ott.optics.beam.Beam`
      %   - xyz (3x1 numeric) -- Dipole position in beam

      % Defer all argument checks to ScatteredBeam
      sbeam = ott.optics.dipole.ScatteredBeam(beam, dipole, xyz);
    end

    function f = force(dipole, beam, xyz)
      % Calculate the total force on a dipole in a beam
      %
      % Usage
      %   f = dipole.force(beam, xyz)
      %
      % Parameters
      %   - beam -- A :class:`ott.optics.beam.Beam`
      %   - xyz (3xN numeric) -- Dipole positions in beam

      % Check inputs
      assert(isa(beam, 'ott.optics.beam.Beam'), ...
        'beam must be a ott.optics.beam.Beam objects');
      assert(isnumeric(xyz) && all(size(xyz, 1) == 3), ...
        'position must be 3xN numeric');

      f = dipole.force_scat(beam, xyz) + dipole.force_grad(beam, xyz);
    end

    function f = force_scat(dipole, beam, xyz)
      % Calculate the scattered force on a dipole in a beam
      %
      % Usage
      %   f = dipole.force_scat(beam, xyz)
      %
      % Parameters
      %   - beam -- A :class:`ott.optics.beam.Beam`
      %   - xyz (3xN numeric) -- Dipole positions in beam

      % Check inputs
      assert(isa(beam, 'ott.optics.beam.Beam'), ...
        'beam must be a ott.optics.beam.Beam objects');
      assert(isnumeric(xyz) && all(size(xyz, 1) == 3), ...
        'position must be 3xN numeric');

      % TODO: Calculate Cpr
      error('Not yet implemented');

      % Calculate Poynting vector
      % TODO: For paraxial beams this should have an efield default
      % and an explicit z dependence
      S = beam.poynting(xyz);

      % Calculate force (Harada Eq 10)
      f = 1./beam.speed .* Cpr .* 0.5.*real(S);
    end

    function f = force_grad(dipole, beam, xyz, varargin)
      % Calculate the gradient force on a dipole in a beam
      %
      % Usage
      %   f = dipole.force_grad(beam, xyz)
      %
      % Parameters
      %   - beam -- A :class:`ott.optics.beam.Beam`
      %   - xyz (3xN numeric) -- Dipole positions in beam
      %
      % Optional named arguments
      %   - dx (numeric) -- step size to use for numerical derivatives.
      %     Default: ``beam.wavelength*1.0e-3``.

      % Check inputs
      assert(isa(beam, 'ott.optics.beam.Beam'), ...
        'beam must be a ott.optics.beam.Beam objects');
      assert(isnumeric(xyz) && all(size(xyz, 1) == 3), ...
        'position must be 3xN numeric');

      ip = inputParser;
      ip.addParameter('dx', 1.0e-3*beam.wavelength);
      ip.parse(varargin{:});

      % Calculate electric field at positions
      E = beam.efield(xyz);

      % Calculate electric field around positions (using forward difference)
      % TODO: We should defer to the beam to calculate this,
      %   this could have a analytical values for certain beams or a
      %   smart thing with a grid of xyz points
      dx = p.Results.dx;
      Ex = beam.efield(xyz + [dx;0;0]);
      Ey = beam.efield(xyz + [0;dx;0]);
      Ez = beam.efield(xyz + [0;0;dx]);
      Eforward = [Ex; Ey; Ez];

      % Calculate polarisation (Harada Eq 9)
      % TODO: Should this be the total or partial polarizability
      if isvector(dipole.polarizability) || isscalar(dipole.polarizability)
        p = dipole.polarizability .* E;
      else
        % Polarizability is 3x3 matrix
        p = dipole.polarizability * E;
      end

      % Calculate force (Harada Eq 13)
      f = p .* (Eforward - E) ./ dx;

    end
  end
end

