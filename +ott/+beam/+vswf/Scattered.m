classdef Scattered < ott.beam.vswf.Bsc
% A Bsc instance describing a scattered beam.
% Inherits from :class:`Bsc`.
%
% An instance of this class is created by the Bsc object when
% a particle scatters a beam.  By default, this object keeps track of the
% incident beam and the T-matrix, allowing easy total-field and
% scattered-field calculation.
%
% Scattered beams can either have a total-field type or scattered-field
% type, depending on if the beam shape coefficients describe only the
% scattered fields or the total fields.
%
% Properties
%   - incident_beam     -- The incident beam that was scattered (or [])
%   - tmatrix           -- The T-matrix which scattered the beam (or [])
%   - type              -- Type of beam (scattered or total)
%
% Methods
%   - totalField        -- Calculate the total field representation
%   - scatteredField    -- Calculated the scattered field representation

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

% TODO: This should be re-written to inherit from ott.beam.Scattered

  properties (Hidden)
    typeInternal       % Value of type parameter
  end

  properties
    incident_beam      % The incident beam that was scattered (or [])
    tmatrix            % The T-matrix which scattered the beam (or [])
  end

  properties (Dependent)
    type               % Type of beam (scattered or total)
  end

  methods
    function beam = Scattered(varargin)
      % Calculate the scattered beam from a beam and T-matrix
      %
      % Usage
      %   beam = Scattered(incident_beam, tmatrix, ...)
      %   Calculate the scattered beam using a T-matrix.
      %
      %   beam = Scattered(a, b, type, ...)
      %   Construct a new Scattered beam with the specified beam shape
      %   coefficients and scattered beam type.
      %
      %   beam = Scattered(bsc, type, ...)
      %   Convert an existing Bsc beam to a Scattered beam.
      %   Does not change the BSC values, just changes the class.
      %
      % Parameters
      %   - incident_beam (Bsc) -- Incident beam object
      %   - tmatrix (Tmatrix) -- T-matrix describing scattering
      %   - type (enum) -- Type of scattered beam.
      %   - a,b (numeric) -- Beam shape vectors.
      %   - bsc (vswf.bsc.Bsc) -- Existing beam to convert.
      %
      % Optional named arguments
      %   - store_tmatrix (logical) -- If false, the T-matrix is not stored.
      %     Default: ``true``.
      %
      %   - store_beam (logical) -- If false, the incident beam is not stored.
      %     Default: ``true``.
      %
      %   - position (3xN numeric) -- Translation applied to the beam
      %     before the beam is scattered by the particle.  Default: ``[]``.
      %     Only used with T-matrix input.
      %
      %   - rotation (3x3N numeric) -- Rotation applied to the beam,
      %     calculates the scattered beam and applies the inverse rotation,
      %     effectively rotating the particle.  Default: ``[]``.
      %     Only used with T-matrix input.
      %
      %   - combine (enum|empty) -- Beam combination method.  Can be
      %     'incoherent' (ignored), 'coherent', or empty (ignored).
      %     Default: ``[]``.  Only used with T-matrix input.
      %
      % If both position and rotation are present, the translation is
      % applied first, followed by the rotation.
      % If both position and rotation are arrays, they must have the same
      % number of locations.

      % Setup input parser
      p = inputParser;
      p.addParameter('store_beam', true);
      p.addParameter('store_tmatrix', true);
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('combine', []);

      % Set-up arguments to pass to base class constructor
      if nargin == 0

        basearg = {};
        type = 'scattered';

        % Parse remaining args
        p.parse();

        % No incident beam or T-matrix
        tmatrix = [];
        incident_beam = [];

      elseif nargin >= 2 && isa(varargin{1}, 'ott.beam.vswf.Bsc') && ...
          isa(varargin{2}, 'ott.optics.vswf.tmatrix.Tmatrix')

        % Parse arguments
        incident_beam = varargin{1};
        tmatrix = varargin{2};
        p.parse(varargin{3:end});

        % Copy the Bsc properties of the incident beam
        basearg = {incident_beam};
        type = [];

      elseif nargin >= 2 && isa(varargin{1}, 'ott.beam.vswf.Bsc')

        % Convert Bsc to Scattered
        basearg = varargin(1);
        type = varargin{2};

        % Parse remaining args
        p.parse(varargin{3:end});

        % No incident beam or T-matrix
        tmatrix = [];
        incident_beam = [];

      elseif nargin >= 3

        % Construct new Bsc for coefficients
        basearg = varargin(1:2);
        type = varargin{3};

        % Parse remaining args
        p.parse(varargin{4:end});

        % No incident beam or T-matrix
        tmatrix = [];
        incident_beam = [];

      else
        error('Invalid number of input arguments');
      end

      % Call base constructor
      beam = beam@ott.beam.vswf.Bsc(basearg{:});

      % If type is set, we are done, otherwise other calculations to do
      if isempty(type)
        % Calculate scattered beam (assigns type from T-matrix)
        beam = beam.calculateScatteredBeam(tmatrix, ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation, ...
            'combine', p.Results.combine);
      else
        beam.typeInternal = type;
      end

      % Store T-matrix (if present/requested)
      if p.Results.store_tmatrix
        beam.tmatrix = tmatrix;
      else
        beam.tmatrix = [];
      end

      % Store beam (if present/requested)
      if p.Results.store_beam
        beam.incident_beam = incident_beam;
      else
        beam.incident_beam = [];
      end
    end

    function beam = totalField(beam, ibeam)
      % Calculate the total field representation of the beam
      %
      % Usage
      %   total_beam = beam.totalField(incident_beam)
      %
      % Parameters
      %   - incident_beam (Bsc) -- Optional if incident_beam is
      %     set internally.

      if nargin == 1
        assert(~isempty(beam.incident_beam), ...
            'incident_beam must be specified in function or class');
        ibeam = beam.incident_beam;
      end

      assert(isa(ibeam, 'ott.beam.vswf.Bsc'), ...
          'incident_beam must be an vswf.bsc.Bsc object');

      switch beam.typeInternal
        case 'total'
          % Nothing to do

        case 'scattered'
          beam = 2*beam + ibeam;
          beam.typeInternal = 'total';

        otherwise
          error('Unsupported beam type');
      end
    end

    function beam = scatteredField(beam, ibeam)
      % Calculate the scattered field representation of the beam
      %
      % Usage
      %   scattered_beam = beam.totalField(incident_beam)
      %
      % Parameters
      %   - incident_beam (Bsc) -- Optional if incident_beam is
      %     set internally.

      if nargin == 1
        assert(~isempty(beam.incident_beam), ...
            'incident_beam must be specified in function or class');
        ibeam = beam.incident_beam;
      end

      assert(isa(ibeam, 'ott.beam.vswf.Bsc'), ...
          'incident_beam must be an vswf.bsc.Bsc object');

      switch beam.typeInternal
        case 'total'
          beam = 0.5*(beam - ibeam);
          beam.typeInternal = 'scattered';

        case 'scattered'
          % Nothing to do

        otherwise
          error('Unsupported beam type');
      end
    end

    function [sbeam, tbeam] = calculateScatteredBeam(tbeam, ...
        tmatrix, varargin)
      % Calculate the scattered beam
      %
      % Usage
      %   [sbeam, tbeam] = beam.calculateScatteredBeam(tmatrix, ...)
      %   Calcualtes the scattered beam ``sbeam`` and the translated
      %   and possibly rotated beam ``tbeam``.
      %
      % For optional arguments, see Bsc.scatter or Scattered constructor.

      p = inputParser;
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('combine', []);
      p.parse(varargin{:});

      % Apply combine parameter
      assert(any(strcmpi(p.Results.combine, {'coherent', 'incoherent'})) ...
          || isempty(p.Results.combine), ...
          'combine must be empty, ''cohernet'' or ''incoherent''');
      if strcmpi(p.Results.combine, 'coherent')
        tbeam = sum(tbeam);
      end

      % Check if we need to calculate multiple scatters
      if size(p.Results.position, 2) > 1 || size(p.Results.rotation, 2) > 3
        [sbeam, tbeam] = ott.utils.prxfun(...
            @(varargin) tbeam.scatter(tmatrix, varargin{:}), 1, ...
            'position', p.Results.position, ...
            'rotation', p.Results.rotation, ...
            'zeros', @(x) repmat(tbeam, x));

        % Join outputs (convert from beam arrays to combined beams)
        sbeam = sbeam(1).append(sbeam(2:end));
        tbeam = tbeam(1).append(tbeam(2:end));

        return;
      end

      % Determine the maximum tmatrix.Nmax(2) and check type
      maxNmax1 = 0;
      maxNmax2 = 0;
      tType = tmatrix(1).type;
      for ii = 1:numel(tmatrix)
        maxNmax1 = max(maxNmax1, tmatrix(ii).Nmax(1));
        maxNmax2 = max(maxNmax2, tmatrix(ii).Nmax(2));
        if ~strcmpi(tmatrix(ii).type, tType)
          error('T-matrices must be same type');
        end
      end

      % If the T is scattered, we can save time by throwing away columns
      if strcmpi(tmatrix(1).type, 'scattered')
        maxNmax2 = min(maxNmax2, tbeam.Nmax);
      end

      % Ensure all T-matrices are the same size
      for ii = 1:numel(tmatrix)
        tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
      end

      % Apply translation to the beam
      if ~isempty(p.Results.position)

        % Requires scattered beam, convert if needed
        if ~strcmpi(tmatrix(1).type, 'scattered')
          maxNmax2 = min(maxNmax2, tbeam.Nmax);
          for ii = 1:numel(tmatrix)
            tmatrix(ii).type = 'scattered';
            tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
          end
        end

        % Apply translation
        % We need Nmax+1 terms for the force calculation
        tbeam = tbeam.translateXyz(p.Results.position, 'Nmax', maxNmax2+1);
      end

      % Apply rotation to the beam
      rbeam = tbeam;
      if ~isempty(p.Results.rotation)
        [rbeam, D] = rbeam.rotate(p.Results.rotation, ...
            'Nmax', maxNmax1);
      end

      % Ensure the Nmax for the inner dimension matches
      if strcmpi(tmatrix(1).type, 'scattered')
        % T-matrix is already done
        rbeam = rbeam.set_Nmax(maxNmax2, 'powerloss', 'ignore');
      else
        for ii = 1:numel(tmatrix)
          tmatrix(ii) = tmatrix(ii).set_Nmax([maxNmax1, rbeam.Nmax], ...
              'powerloss', 'ignore');
        end
        if ~strcmpi(tmatrix(1).type, 'internal')
          ott.warning('ott:Bsc:scatter', ...
              'It may be more optimal to use a scattered T-matrix');
        end
      end

      % Calculate the resulting beams
      sbeam = ott.beam.vswf.Bsc();
      for ii = 1:numel(tmatrix)
        sbeam = [sbeam, tmatrix(ii).data * rbeam];
      end

      % Apply the inverse rotation
      if ~isempty(p.Results.rotation)

        % This seems to take a long time
        %sbeam = sbeam.rotate('wigner', D');

        sbeam = sbeam.rotate(inv(p.Results.rotation));
      end

      % Assign a type to the resulting beam
      switch tmatrix(1).type
        case 'total'
          sbeam.typeInternal = 'total';
          sbeam.basis = 'regular';
        case 'scattered'
          sbeam.typeInternal = 'scattered';
          sbeam.basis = 'outgoing';
        case 'internal'
          sbeam.typeInternal = 'internal';
          sbeam.basis = 'regular';

          % Wavelength has changed, update it
          sbeam.wavenumber = tmatrix(1).wavenumber_particle;

        otherwise
          error('Unrecognized T-matrix type');
      end
    end
  end

  methods
    function beam = set.typeInternal(beam, val)
      % Set the internal type parameter

      assert(any(strcmpi(val, {'scattered', 'total', 'internal'})), ...
          'type must be one of ''scattered'' or ''total'' or ''internal''');

      beam.typeInternal = val;
    end

    function val = get.type(beam)
      % Get the internal type value
      val = beam.typeInternal;
    end
    function beam = set.type(beam, val)
      % Set the beam type, checking it is a valid type first

      assert(any(strcmpi(val, {'scattered', 'total', 'internal'})), ...
          'type must be one of ''scattered'' or ''total'' or ''internal''');

      assert(~isempty(beam.incident_beam), ...
          'Need incident beam to change type');

      % Change the beam
      if strcmpi(val, 'total')
        beam = beam.totalField();
      elseif strcmpi(val, 'scattered')
        beam = beam.scatteredField();
      else
        error('Cannot convert to internal beam from external beam');
      end
    end

    function beam = set.tmatrix(beam, val)
      assert(isempty(val) || isa(val, 'ott.optics.vswf.tmatrix.Tmatrix'), ...
          'tmatrix must be a valid vswf.tmatrix.Tmatrix');
      beam.tmatrix = val;
    end

    function beam = set.incident_beam(beam, val)
      assert(isempty(val) || isa(val, 'ott.beam.vswf.Bsc'), ...
          'incident_beam must be a vswf.bsc.Bsc');
      beam.incident_beam = val;
    end
  end
end
