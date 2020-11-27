classdef BscBeam < ott.beam.ArrayType & ott.beam.properties.IndexOmegaProps
% Beam class encapsulating a BSC instance.
% Inherits from :class:`+ott.+beam.ArrayType`.
%
% This class stores an internal :class:`+ott.+bsc.Bsc` instance which it
% uses for calculating fields.
%
% Methods in this class assume the beam is a regular beam (i.e., not
% the outgoing fields from a scattered beam).
%
% Properties
%   - data          -- Internal BSC instance describing beam
%   - apparentPower -- Apparent power of the BSC data
%
% Methods
%   - recalculate   -- Can be overloaded by sub-classes to update data
%   - efield, hfield, ... -- Calculate fields in SI units, uses the
%     :class:`ott.bsc.Bsc` methods internally.
%   - force, torque, spin -- Calculate force/torque/spin in SI units
%
% Supported casts
%   - ott.bsc.Bsc   -- Get the BSC data after applying transformations
%
% Additional properties/methods inherited from :class:`Beam`.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    data      % Internal BSC instance describing beam
  end

  properties (Dependent)
    apparentPower   % Apparent power of the BSC data
  end
  
  methods (Static)
    function DefaultScatterProgressCallback(data)
      % Default progress callback for :meth:`scatterBsc`.
      %
      % Prints the progress to the terminal.
      %
      % Usage
      %   DefaultScatterProgressCallback(data)
      %
      % Parameters
      %   - data (struct) -- Structure with two fields: ``iteration``
      %     and ``total``.
      
      % Print progress every 5%
      iout = round(0.1*data.total);
      if mod(data.iteration, iout) == 0
        disp(['Scatter progress... ' ...
          num2str((data.iteration)/data.total*100) '%']);
      end
    end
  end

  methods
    function bm = BscBeam(varargin)
      % Construct a Bsc beam instance.
      %
      % Usage
      %   beam = BscBeam(data, ...)
      %
      % Optional named parameters
      %   - data (ott.bsc.Bsc) -- Initial data for beam or empty.
      %     Default: ``ott.bsc.Bsc.empty()``.
      %
      %   - arrayType (enum) -- The array type for Bsc beams containing
      %     multiple bsc data.  Default: ``'coherent'``.
      %
      %   - index_medium (numeric) -- Refractive index of the medium.
      %     Default: ``1.0``.
      %
      %   - omega (numeric) -- Optical angular frequency [Hz].
      %     Default: ``3e8/1064e-9*2*pi`` (i.e., default vacuum
      %     wavelength is 1064 nm).

      p = inputParser;
      p.addOptional('data', ott.bsc.Bsc.empty(), @(x) isa(x, 'ott.bsc.Bsc'));
      p.addParameter('arrayType', 'coherent');
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      [omega, index, unmatched] = ott.beam.properties. ...
          IndexOmegaProps.parseArgs(unmatched{:});

      bm = bm@ott.beam.properties.IndexOmegaProps(omega, index);
      bm = bm@ott.beam.ArrayType(...
        'arrayType', p.Results.arrayType, unmatched{:});
      bm.data = p.Results.data;
    end

    function bsc = ott.bsc.Bsc(beam, Nmax)
      % Get the BSC data for a specific Nmax.
      %
      % If the current data is empty or the current Nmax is less than
      % the desired Nmax, calls :meth:`recalculate` and then returns the
      % internal data.
      %
      % If position and rotation have been set, rotations and translations
      % are applied to the beam (rotations first, translations second).
      % Translations of this beam type require a larger internal Nmax,
      % if the Nmax is not large enough for the translation,
      % :meth:`recalculate` is called.
      %
      % Usage
      %   bsc = ott.bsc.Bsc(beam, Nmax)
      %
      % Parameters
      %   - beam (ott.beam.BscBeam) -- The beam to get data from.
      %
      %   - Nmax (numeric) -- Desired Nmax.  The default is either 0
      %     or the current beam Nmax minus the Nmax required for any
      %     translations, whichever is greater.

      % Calculate Nmax required after translation
      % This assumes the beam is a regular beam
      % Also assumes generic translation
      dist = vecnorm(beam.position);
      tNmax = ott.utils.ka2nmax(beam.wavenumber*dist);

      if nargin < 2
        Nmax = max(0, max([0, beam.data.Nmax]) - tNmax);
      end
      assert(isnumeric(Nmax) && isscalar(Nmax), ...
        'Nmax must be numeric scalar');

      if isempty(beam.data) || max([0, beam.data.Nmax]) < (Nmax + tNmax)
        if Nmax + tNmax == 0
          beam.data = ott.bsc.Bsc();
        else
          beam = beam.recalculate(Nmax + tNmax);
        end
      end

      bsc = beam.data;

      % Rotate the beam
      bsc = bsc.rotate(beam.rotation);

      % Translate the beam
      % Assumes beam is a regular beam
      bsc = bsc.translateXyz(-beam.position ./ beam.wavelength, 'Nmax', Nmax);

      % Combine beams now if coherent
      if strcmpi(beam.arrayType, 'coherent')
        bsc = sum(subsum(bsc));
      end

      % Apply scale to data
      bsc = bsc * beam.scale;
    end

    function beam = recalculate(beam, ~) %#ok<INUSD>
      % Re-calculate data for a specified Nmax
      %
      % This implementation simply raises an error when called.
      % This function should be overloaded by sub-classes of BscBeam
      % to implement the desired behaviour.
      %
      % Usage
      %   beam = beam.recalcualte(Nmax)
      %
      % Parameters
      %   - Nmax (numeric) -- Desired Nmax for beam.

      error('ott:beam:BscBeam:recalculate_not_implemented', ...
          'recalculate is not implemented for these type of beam');
    end

    function sbeam = scatterInternal(ibeam, particle)
      % Calculate how a particle scatters the beam
      %
      % Usage
      %   sbeam = scatter(ibeam, particle)
      %
      % Returns
      %   - sbeam (ott.beam.Scattered) -- Scattered beam encapsulating
      %     the particle, incident beam, scattered beams(s).  For a
      %     method which doesn't create a scattered beam, see
      %     :meth:`scatterBsc`.
      %
      % Parameters
      %   - particle (ott.particle.Particle) -- Particle with
      %     T-matrix properties (possibly internal and external).
      
      % Work in the particle reference frame
      scat_position = particle.position;
      scat_rotation = particle.rotation;
      particle.position = [0;0;0];
      particle.rotation = eye(3);
      ibeam = ibeam.translateXyz(-scat_position).rotate(scat_rotation.');
      
      % Calculate scattered Bscs
      [ibsc, sbsc] = ibeam.scatterBsc(particle);
      
      % Calculate external component
      if ~isempty(sbsc)
        sbsc = sbsc.rotate(scat_rotation);
        sbeam = ott.beam.BscOutgoing(sbsc, ...
          'index_medium', ibeam.index_medium, 'omega', ibeam.omega);
        sbeam.rotation = scat_rotation.';
      else
        sbeam = [];
      end
      
      % Calculate internal component
      if ~isempty(particle.tinternal) && isa(particle.tinternal, 'ott.tmatrix.Homogeneous')
        sint = particle.tinternal * ibsc;
        sint = sint.rotate(scat_rotation);
        index_particle = ibeam.index_medium * particle.tinternal.index_relative;
        nbeam = ott.beam.BscBeam(sint, ...
          'index_medium', index_particle, 'omega', ibeam.omega);
        nbeam.rotation = scat_rotation.';
      else
        nbeam = [];
      end

      % Package output
      sbeam = ott.beam.Scattered(...
          'scattered', sbeam, 'incident', ibeam, ...
          'particle', particle, 'internal', nbeam, ...
          'position', scat_position, 'rotation', scat_rotation);
    end
    
    function [ibsc, sbsc] = scatterBsc(beam, particle, varargin)
      % Calculate the beam shape coefficients from scattering
      %
      % Usage
      %   [ibsc, sbsc] = beam.scatterBsc(particle, ...)
      %
      % Returns
      %   - ibsc (ott.bsc.Bsc) -- Incident beam shape coefficients
      %     after applying rotations and translations.
      %
      %   - sbsc (ott.bsc.Bsc) -- Scattered beam shape coefficients
      %     in particle reference frame.
      %
      % Parameters
      %   - particle (ott.particle.Particle) -- Scattering particle.
      %     Uses the position, rotation and T-matrix properties.
      %
      % Optional named arguments
      %   - position (3xN numeric | 3x1 cell) -- Particle position.
      %     Can either be Cartesian position vectors or a cell array
      %     with ``{X, Y, Z}`` cartesian coordinates.  Each matrix in
      %     cell array must have the same size.
      %     Default: ``particle.position``.
      %
      %   - rotation (3x3N numeric) -- Particle rotation.
      %     Default: ``particle.rotation``.
      %
      %   - progress (function_handle) -- Function to call for progress
      %     updates during method evaluation.  Takes one argument, see
      %     :meth:`DefaultScatterProgressCallback` for more information.
      %     Default: ``[]`` (for N < 200) and
      %     ``@DefaultScatterProgressCallback`` (otherwise).
      
      p = inputParser;
      p.addParameter('position', [particle.position]);
      p.addParameter('rotation', [particle.rotation]);
      p.addParameter('progress', []);
      p.parse(varargin{:});
      
      % Get position and rotation
      position = p.Results.position;
      rotation = p.Results.rotation;
      
      % Convert cell array to vector format
      if iscell(position)
        position = [position{1}(:), position{2}(:), position{3}(:)].';
      end
      
      % Check sizes
      assert(isnumeric(position) && ismatrix(position) && size(position, 1) == 3, ...
        'position must be 3xN numeric matrix');
      assert(isnumeric(rotation) && ismatrix(rotation) ...
        && size(rotation, 1) == 3 && mod(size(rotation, 2), 3) == 0, ...
        'rotation must be 3x3N numeric matrix');
      
      % Count amount of work, repmat position/rotation if required
      Npos = size(position, 2);
      Nrot = size(rotation, 2)/3;
      Npart = numel(particle);
      Nwork = max([Npos, Nrot, Npart]);
      assert((Npos == Nwork || Npos == 1) ...
        && (Nrot == Nwork || Nrot == 1) ...
        && (Npart == Nwork || Npart == 1), ...
        'Number of position/rotation/particles must be equal or 1');
      if Npos ~= Nwork, position = repmat(position, [1, Nwork]); end
      if Nrot ~= Nwork, rotation = repmat(rotation, [1, Nwork]); end
      if Npart ~= Nwork, particle = repmat(particle, [1, Nwork]); end

      % Handle default argument for progress callback
      progress_cb = p.Results.progress;
      if isempty(progress_cb)
        if Nwork >= 200
          progress_cb = @ott.beam.BscBeam.DefaultScatterProgressCallback;
        else
          progress_cb = @(x) [];
        end
      end

      % Get required Nmax for beam data
      Nmax = 0;
      if numel([particle.tmatrix]) == numel(particle)
        ts = [particle.tmatrix];
        Nmax = max(cat(1, [0, 0], ts.Nmax), 1);
        Nmax = Nmax(2);
      end
      
      % Pre-allocate space for results
      ibsc(Nwork) = ott.bsc.Bsc();
      if numel([particle.tmatrix]) == numel(particle)
        sbsc(Nwork) = ott.bsc.Bsc();
      else
        sbsc = [];
      end
        
      % Report initial progress
      progress_cb(struct('total', Nwork, 'iteration', 0));
      
      for ii = 1:Nwork

        % Apply particle position to beam and particle
        % Make a copy of the beam data (new beam uses original incidnet beam)
        tbeam = beam.translateXyz(-position(:, ii));

        % Get bsc data
        ibsc(ii) = ott.bsc.Bsc(tbeam, Nmax);

        % Apply particle rotation
        % We are now in the particle frame with both ibsc and sbsc
        ibsc(ii) = ibsc(ii).rotate(rotation(:, (1:3) + (ii-1)*3).');

        % Calculate external component
        if ~isempty(particle(ii).tmatrix)
          sbsc(ii) = particle(ii).tmatrix * ibsc(ii);
        end
        
        % Report progress
        progress_cb(struct('total', Nwork, 'iteration', ii));
      end
      
      if iscell(p.Results.position)
        ibsc = reshape(ibsc, size(p.Results.position{1}));
        sbsc = reshape(sbsc, size(p.Results.position{1}));
      end
    end

    %
    % Field calculation functions
    %

    function varargout = efarfield(beam, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, H, ...] = beam.efarfield(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      bsc = ott.bsc.Bsc(beam);
      [varargout{1:nargout}] = bsc.efarfield(varargin{:});

      if nargout >= 1
        varargout{1} = varargout{1} .* sqrt(beam.impedance);
      end
    end

    function [H, vswfData] = hfarfield(beam, rtp, varargin)
      % Calculate H far-field (uses :meth:`efarfield`)
      %
      % Usage
      %   [H, data] = beam.hfarfield(rtp, ...)
      %
      % See :meth:`efarfield` for further details.

      [~, H, vswfData] = beam.ehfarfield(rtp, varargin{:});
    end

    function [E, H, data] = ehfarfield(beam, rtp, varargin)
      % Calculate E and H far-fields
      %
      % Usage
      %   [E, H, data] = beam.ehfarfield(rtp, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - rtp (3xN | 2xN numeric) -- Spherical coordinates for field
      %     calculation. Packaged [r; theta; phi] or [theta; phi].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      % Ensure rtp size is 3xN
      [~, rtp] = ott.utils.rtpFarfield(rtp);

      % Calculate E as normal
      [E, data] = beam.efarfield(rtp, varargin{:});

      % Swap theta and phi components
      H = ott.utils.FieldVectorSph(...
          -1i .* E.vrtp([1, 3, 2], :), rtp);
      H = H ./ beam.impedance;
    end

    function varargout = efield(beam, xyz, varargin)
      % Calculate E field from Cartesian coordinates.
      %
      % Usage
      %   [E, data] = beam.efield(xyz, ...)
      %   E is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.efieldRtp(rtp, varargin{:});
    end

    function varargout = hfield(beam, xyz, varargin)
      % Calculate H field from Cartesian coordinates.
      %
      % Usage
      %   [H, data] = beam.hfield(xyz, ...)
      %   H is of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.hfieldRtp(rtp, varargin{:});
    end

    function varargout = ehfield(beam, xyz, varargin)
      % Calculate E and H fields from Cartesian coordinates.
      %
      % Usage
      %   [E, H, data] = beam.ehfield(xyz, ...)
      %   E and H are of type :class:`ott.utils.FieldVector`.
      %
      % Parameters
      %   - xyz (3xN numeric) -- Cartesian coordinates for field calculation.
      %     Units of wavelength.  Packaged [x; y; z].
      %
      % Optional named parameters
      %   - data (ott.utils.VswfData) -- Field data for repeated field
      %     calculation.  Default is an empty VswfData structure.

      % This overloads the default in order to avoid two xyz2rtp casts

      rtp = ott.utils.xyz2rtp(xyz);
      [varargout{1:nargout}] = beam.ehfieldRtp(rtp, varargin{:});
    end

    function varargout = efieldRtp(beam, rtp, varargin)
      % Calculate the electric field (in SI units)
      %
      % Usage
      %   [E, ...] = beam.efieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      % Translate coordinates instead of beam
      rtp = beam.global2localRtp(rtp);
      beam.position = [0;0;0];
      beam.rotation = eye(3);

      bsc = ott.bsc.Bsc(beam);
      [varargout{1:nargout}] = bsc.efieldRtp(...
          rtp./[beam.wavelength;1;1], varargin{:});

      if nargout >= 1
        varargout{1} = varargout{1} ./ sqrt(beam.index_medium) ...
          .* 2 .* beam.wavenumber .* sqrt(beam.impedance);
      end
    end

    function varargout = hfieldRtp(beam, rtp, varargin)
      % Calculate the magnetic field (in SI units)
      %
      % Usage
      %   [H, ...] = beam.hfieldRtp(rtp, ...)
      %
      % Parameters
      %   - rtp (3xN numeric) -- Spherical coordinates for field calculation.
      %     Units of meters/radians.  Packaged [r;theta;phi].
      %
      % See :class:`+ott.+bsc.Bsc` for further details.

      % Translate coordinates instead of beam
      rtp = beam.global2localRtp(rtp);
      beam.position = [0;0;0];
      beam.rotation = eye(3);

      bsc = ott.bsc.Bsc(beam);
      [varargout{1:nargout}] = bsc.hfieldRtp(...
          rtp./[beam.wavelength;1;1], varargin{:});
        
      if nargout >= 1
        varargout{1} = varargout{1} ./ sqrt(beam.index_medium) ...
          .* 2 .* beam.wavenumber ./ sqrt(beam.impedance);
      end
    end

    function F = force(ibeam, sbeam, varargin)
      % Calculate force (in Newtons)
      %
      % Usage
      %   force = incident_beam.force(scattered_beam) -- Calculates
      %   difference in momentum between two beams: F = Scat - Inc.
      %
      %   force = incident_beam.force(particle, ...) -- First calculates
      %   scattering, then calculates force from difference in momentum.
      %   Additional parameters are passed to :meth:`scatterBsc`.
      
      % Calculate force
      F = forceHelper(ibeam, sbeam, @(i, s) i.force(s), varargin{:});

      % Convert units of output to SI
      F = F ./ ibeam.speed;
    end

    function T = torque(ibeam, sbeam, varargin)
      % Calculate torque (in Newton meters)
      %
      % Usage
      %   torque = incident_beam.torque(scattered_beam) -- Calculates
      %   difference in momentum between two beams.
      %
      %   torque = incident_beam.torque(particle, ...) -- First calculates
      %   scattering, then calculates torque from difference in momentum.
      %   Additional parameters are passed to :meth:`scatter`.

      % Calculate force
      T = forceHelper(ibeam, sbeam, @(i, s) i.torque(s), varargin{:});

      % Convert units of output to SI
      T = T ./ ibeam.omega;
    end

    function T = spin(ibeam, sbeam, varargin)
      % Calculate spin (in Newton meters)
      %
      % Usage
      %   spin = incident_beam.spin(scattered_beam) -- Calculates
      %   difference in momentum between two beams.
      %
      %   spin = incident_beam.spin(particle, ...) -- First calculates
      %   scattering, then calculates spin from difference in momentum.
      %   Additional parameters are passed to :meth:`scatter`.

      % Calculate force
      T = forceHelper(ibeam, sbeam, @(i, s) i.spin(s), varargin{:});

      % Convert units of output to SI
      T = T ./ ibeam.omega;
    end
  end

  methods (Hidden)
    function O = forceHelper(ibeam, sbeam, func, varargin)
      % Helper for force/torque methdods
      
      % Calculate scattering if required
      if ~isa(sbeam, 'ott.beam.Beam')
      
        p = inputParser;
        p.addParameter('rotation', [sbeam.rotation]);
        p.KeepUnmatched = true;
        p.parse(varargin{:});
        
        % Calculate scattered field (incident-scattered)
        [ibsc, sbsc] = ibeam.scatterBsc(sbeam, varargin{:});
        
        % Get particle rotation
        particle_rot = p.Results.rotation;
        
        % Convert to incoming-outgoing
        sbsc = ibsc + 2*sbsc;
      else
        
        % Get particle rotation (to be applied later)
        particle_rot = sbeam.rotation;
        
        % Apply the scattered beams translation to ourselves
        % Perhaps this isn't the best thing to do, but I'm not sure
        % where else we should do this.
        ibeam = ibeam.translateXyz(-sbeam.position).rotate(sbeam.rotation.');
        sbeam.position = [0;0;0];
        sbeam.rotation = eye(3);
        
        sbsc = ott.bsc.Bsc(sbeam);
        
        % Use Nmax from sbsc, we don't apply translations to sbsc so
        % this Nmax should correspond to the scattered Nmax.
        ibsc = ott.bsc.Bsc(ibeam, [sbsc.Nmax]);
      end

      % Calculate force using internal methods
      O = func(ibsc, sbsc);
      O = reshape(O, 3, [], numel(ibsc));
      
      % Apply particle rotation term to result
      % This avoid calculating a wigner matrix
      if size(particle_rot, 2) == 3
        O = reshape(particle_rot * O(:, :), size(O));
      else
        for ii = 1:size(particle_rot, 2)/3
          O(:, :, ii) = particle_rot(:, (1:3) + (ii-1)*3) * O(:, :, ii);
        end
      end
    end
    
    function val = defaultVisRangeInternal(beam)
      val = [1,1] .* ott.utils.nmax2ka(...
          max([1, beam.data.Nmax]))./beam.wavenumber ./ sqrt(2);
    end
  end

  methods % Getters/setters
    function beam = set.data(beam, val)
      if isempty(val)
        val = ott.bsc.Bsc.empty();
      end
      assert(isa(val, 'ott.bsc.Bsc'), ...
          'data must be a ott.bsc.Bsc instance or empty');
      beam.data = val;
    end

    function p = get.apparentPower(beam)
      if isempty(beam.data)
        p = nan;
      else
        p = beam.data.power;
      end
    end
  end
end
