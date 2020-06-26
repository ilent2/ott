classdef Bsc < ott.beam.Beam ...
    & ott.beam.properties.AnyArrayType ...
    & ott.beam.properties.VariableMedium
% Class representing beam shape coefficients.
% Inherits from :class:`ott.beam.Beam` 
% and :class:`ott.beam.properties.AnyArrayType`.
%
% Properties
%   - a             -- Beam shape coefficients a vector
%   - b             -- Beam shape coefficients b vector
%   - basis         -- VSWF beam basis (incoming, outgoing or regular)
%   - absdz         -- Absolute cumulative distance the beam has moved
%   - Nmax          -- Truncation number for VSWF coefficients
%
% Inherited properties
%   - power         -- The power of the beam (may be infinite)
%   - omega         -- Beam optical frequency
%   - medium        -- Medium where beam is propagating
%   - position      -- Position of the beam or array
%   - rotation      -- Rotation of the beam or array
%
% Methods
%   - sum         --  Merge the BSCs for the beams contained in this object
%   - translateZ  --  Translates the beam along the z axis
%   - translateXyz -- Translation to xyz using rotations and z translations
%   - translateRtp -- Translation to rtp using rotations and z translations
%   - farfield     -- Calculate fields in farfield
%   - emFieldXyz   -- Calculate field values in cartesian coordinates
%   - emFieldRtp   -- Calculate field values in spherical coordinates
%   - getCoefficients -- Get the beam coefficients [a, b]
%   - getModeIndices -- Get the mode indices [n, m]
%   - visualise      -- Generate a visualisation of the beam near-field
%   - visualiseFarfield -- Generate a visualisation of the beam far-field
%   - visualiseFarfieldSlice  -- Generate scattering slice at specific angle
%   - visualiseFarfieldSphere -- Generate spherical surface visualisation
%   - intensityMoment -- Calculate moment of beam intensity in the far-field
%   - force       --  Calculate change in linear momentum between beams
%   - torque      --  Calculate change in angular momentum between beams
%   - spin        --  Calculate change in spin between beams
%   - size        -- Get the size of the beam array
%
% Static methods:
%   - FromDenseBeamVectors -- Build Bsc from Dense BSC beam vectors
%
% See also Bsc, ott.BscPmGauss, ott.BscPlane.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% TODO: Review documentation
% TODO: Should default basis be regular and automatic conversion for
%   methods that use far-field (unless overridden)

  properties (SetAccess=protected)
    a           % Beam shape coefficients a vector
    b           % Beam shape coefficients b vector
  end

  properties
    absdz       % Absolute cumulative distance the beam has moved
    basis       % VSWF beam basis (incoming, outgoing or regular)
  end

  properties (Dependent)
    Nmax        % Truncation number for VSWF coefficients
    power       % Power of the beam shape coefficients
  end

  methods (Static)
    function bsc = FromDenseBeamVectors(a, b, n, m)
      % Construct a new Bsc instance from dense beam vectors
      %
      % Usage
      %   bsc = Bsc.FromDenseBeamVectors(a, b, n, m)
      %
      % Parameters
      %   - a,b (numeric) -- Dense beam vectors.  Can be Nx1 or NxM for
      %     M beams.  Each row corresponds to n/m indices.
      %
      %   - n,m (numeric) -- Mode incites for beam shape coefficients.
      %     Must have same number of rows as a/b.

      % Check inputs
      assert(all(size(a) == size(b)), 'size of a must match size of b');
      assert(numel(n) == numel(m), 'length of m and n must match');
      assert(numel(n) == size(a, 1), 'length of n must match rows of a');

      % Generate combine indices
      ci = ott.utils.combined_index(n, m);

      % Calculate total_order and number of beams
      nbeams = size(a, 2);
      Nmax = max(n);
      total_orders = ott.utils.combined_index(Nmax, Nmax);

      % Replicate ci for 2-D beam vectors (multiple beams)
      [ci, cinbeams] = ndgrid(ci, 1:nbeams);

      if ~isempty(ci)
        fa = sparse(ci, cinbeams, a, total_orders, nbeams);
        fb = sparse(ci, cinbeams, b, total_orders, nbeams);
      else
        fa = sparse(0, 0);
        fb = sparse(0, 0);
      end

      bsc = ott.beam.vswf.Bsc(fa, fb);
    end

    function args = likeProperties(other, args)
      % Construct an array of like-properties
      if isa(other, 'ott.beam.vswf.Bsc')
        args = ott.utils.addDefaultParameter('a', other.a, args);
        args = ott.utils.addDefaultParameter('b', other.b, args);
        args = ott.utils.addDefaultParameter('basis', other.basis, args);
        args = ott.utils.addDefaultParameter('absdz', other.absdz, args);
      end
      args = ott.beam.Beam.likeProperties(other, args);
      args = ott.beam.properties.AnyArrayType.likeProperties(other, args);
      args = ott.beam.properties.VariableMedium.likeProperties(other, args);
    end

    function beam = like(other, varargin)
      % Construct a VSWF beam like another beam
      args = ott.beam.vswf.Bsc.likeProperties(other, varargin);
      beam = ott.beam.vswf.Bsc(args{:});
    end

    function beam = empty(varargin)
      % Create an empty Bsc object
      %
      % Usage
      %   beam = Bsc.empty()
      %
      %   beam = Bsc.empty(sz, ...)
      %
      % Parameters
      %   - sz (numeric) -- Number of beams in empty array.
      %     Can either be [1, num] or [num].
      %
      % Additional arguments are passed to constructor.

      p = inputParser;
      p.addOptional('sz', 0, @isnumeric);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      sz = p.Results.sz;

      assert(isnumeric(sz) && (sz(1) == 1 || isscalar(sz)), ...
        'sz must be numeric scalar or first element must be 1');

      if numel(sz) == 2
        sz = sz(2);
      end

      a = zeros(0, sz);
      b = zeros(0, sz);

      beam = ott.beam.vswf.Bsc(a, b, unmatched{:});
    end
  end

  methods
    function beam = Bsc(varargin)
      % Construct a new beam object
      %
      % Usage
      %   beam = Bsc() Constructs an empty Bsc beam.
      %
      %   beam = Bsc(a, b, ...) constructs beam from a/b coefficients.
      %
      % Parameters
      %   - a,b (numeric) -- Vectors of VSWF coefficients
      %
      % Optional named arguments
      %   - basis (enum) -- VSWF basis: incoming, outgoing or regular.
      %     Default: ``'regular'``.
      %
      %   - absdz (numeric) -- Initial displacement of the beam.  This is
      %     used to keep track of when the beam may be displaced too far.
      %     Default: ``0``.
      %
      %   - array_type (enum) -- Type of beam array.  Can be one of
      %     'coherent', 'array', or 'incoherent'.  Default: ``'coherent'``.
      %
      % Unmatched arguments are passed to the base class.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('a', [], @isnumeric);
      p.addOptional('b', [], @isnumeric);
      p.addParameter('basis', 'regular');
      p.addParameter('absdz', 0.0);
      p.addParameter('array_type', 'coherent');
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Call base classes
      beam = beam@ott.beam.properties.AnyArrayType(p.Results.array_type);
      beam = beam@ott.beam.Beam(unmatched{:});

      % Store properties
      beam = beam.setCoefficients(p.Results.a, p.Results.b);
      beam.basis = p.Results.basis;
      beam.absdz = p.Results.absdz;
    end

    function beam = ott.beam.vswf.Bsc(varargin)
      % Cast to vswf.FarfieldPm
      beam = ott.beam.vswf.Bsc.like(varargin{:});
    end

    function bsc = ott.beam.vswf.Array(bsc)
      % Split the BSC into an array of BSC

      array_type = bsc.array_type;

      beam_data = cell(1, numel(bsc));
      for ii = 1:numel(bsc)
        beam_data{ii} = bsc(ii);
      end

      bsc = ott.beam.vswf.Array(array_type, beam_data);
    end

    function beam = defaultArrayType(beam, array_type, elements)
      % Construct a new array for this type

      if beam.contains('array')
        assert(strcmpi(array_type, 'array'), ...
            'type must be array for beams containing generic arrays');
      elseif beam.contains('incoherent')
        assert(strcmpi(array_type, {'array', 'incoherent'}), ...
            'type must be array/incoherent for incoherent content');
      end

      if nargin == 2
        beam = @(arg) ott.beam.vswf.Array(array_type, arg);
      else
        beam = ott.beam.vswf.Array(array_type, elements);
      end
    end

    function varargout = size(beam, varargin)
      % Get the number of beams contained in this object
      %
      % Usage
      %   sz = size(beam)   or    sz = beam.size()
      %   For help on arguments, see builtin ``size``.
      %
      % The leading dimension is always 1.  May change in future.

      sz = size(beam.a);
      sz(1) = 1;

      [varargout{1:nargout}] = ott.utils.size_helper(sz, varargin{:});
    end

    function beam = repmat(beam, varargin)

      dims = [varargin{:}];
      assert(dims(1) == 1 && length(dims) == 2, 'Bsc beam array must be 1xN');

      % Repeat the beam coefficients
      beam = beam.setCoefficients(...
        repmat(beam.a, dims), repmat(beam.b, dims));
    end

    function [E, H, data] = ehfarfield(beam, rtp, varargin)
      %FARFIELD finds far field at locations theta, phi.
      %
      % [E, H] = beam.farfield(theta, phi) calculates the farfield
      % at locations theta, phi.  Returns 3xN matrices of the fields
      % in spherical coordinates (r, t, p), the radial component is zero.
      %
      % Theta is the rotation off the z-axis, phi is the rotation about
      % the z-axis.
      %
      % [E, H, data] = beam.farfield(theta, phi, 'saveData', true) outputs
      % a matrix of data the can be used for repeated calculation.
      %
      % Optional named arguments:
      %    'calcE'   bool   calculate E field (default: true)
      %    'calcH'   bool   calculate H field (default: nargout == 2)
      %    'saveData' bool  save data for repeated calculation (default: false)
      %    'data'    data   data saved for repeated calculation.
      %
      % If either calcH or calcE is false, the function still returns
      % E and H as matricies of all zeros.

      ip = inputParser;
      ip.addParameter('calcE', true);
      ip.addParameter('calcH', nargout >= 2);
      ip.addParameter('saveData', false);
      ip.addParameter('data', []);
      ip.parse(varargin{:});

      if size(rtp, 1) == 3
        theta = rtp(2, :).';
        phi = rtp(3, :).';
      else
        theta = rtp(1, :);
        phi = rtp(2, :);
      end

      [theta,phi] = ott.utils.matchsize(theta,phi);

      [theta_new,~,indY]=unique(theta);
      [phi_new,~,indP]=unique(phi);

      Etheta=zeros(length(theta),1);
      Ephi=zeros(length(theta),1);

      Htheta=zeros(length(theta),1);
      Hphi=zeros(length(theta),1);

      if strcmp(beam.basis, 'incoming')

        a = beam.a;
        b = beam.b;
        p = zeros(size(beam.a));
        q = zeros(size(beam.b));

      elseif strcmp(beam.basis, 'outgoing')

        a = zeros(size(beam.a));
        b = zeros(size(beam.a));
        p = beam.a;
        q = beam.b;

      else

        error('Regular wavefunctions go to zero in far-field');

      end

      a = ott.utils.threewide(a);
      b = ott.utils.threewide(b);
      p = ott.utils.threewide(p);
      q = ott.utils.threewide(q);

      [n,m]=ott.utils.combined_index(find(abs(beam.a)|abs(beam.b)));

      % Alocate memory for output data
      data = [];
      if ip.Results.saveData
        data = zeros(numel(indY), 0);
      end

      % Start a counter for accessing the data
      if ~isempty(ip.Results.data)
        dataCount = 0;
      end

      for nn = 1:max(n)

        vv=find(n==nn);
        if isempty(vv)
          continue;
        end

        %this makes the vectors go down in m for n.
        % has no effect if old version code.
        Nn = 1/sqrt(nn*(nn+1));

        % Create index arrays for a, b, q, p
        index=nn*(nn+1)+m(vv);
        aidx = full(a(index));
        bidx = full(b(index));
        pidx = full(p(index));
        qidx = full(q(index));

        if isempty(ip.Results.data)

          [~,Ytheta,Yphi] = ott.utils.spharm(nn,m(vv), ...
              theta_new,zeros(size(theta_new)));

          [PHI,M]=ndgrid(phi_new, m(vv));

          expimphi=exp(1i*M.*PHI);

          % Create full matrices (opt, R2018a)
          YthetaExpf = Ytheta(indY, :).*expimphi(indP, :);
          YphiExpf = Yphi(indY, :).*expimphi(indP, :);

          % Save the data if requested
          if ip.Results.saveData
            data(:, end+(1:size(Ytheta, 2))) = YthetaExpf;
            data(:, end+(1:size(Ytheta, 2))) = YphiExpf;
          end

        else

          % Load the data if present
          YthetaExpf = ip.Results.data(:, dataCount+(1:length(vv)));
          dataCount = dataCount + length(vv);
          YphiExpf = ip.Results.data(:, dataCount+(1:length(vv)));
          dataCount = dataCount + length(vv);

        end

        % Now we use full matrices, we can use matmul (opt, R2018a)
        if ip.Results.calcE
          Etheta = Etheta + Nn * ...
            ( YphiExpf*((1i)^(nn+1)*aidx + (-1i)^(nn+1)*pidx) ...
            + YthetaExpf*((1i)^nn*bidx + (-1i)^nn*qidx) );
          Ephi = Ephi + Nn * ...
            (-YthetaExpf*((1i)^(nn+1)*aidx + (-1i)^(nn+1)*pidx) ...
            + YphiExpf*((1i)^nn*bidx + (-1i)^nn*qidx) );
        end

        if ip.Results.calcH
          Htheta = Etheta + Nn * ...
            ( YphiExpf*((1i)^(nn+1)*bidx + (-1i)^(nn+1)*qidx) ...
            + YthetaExpf*((1i)^nn*aidx + (-1i)^nn*pidx) );
          Hphi = Ephi + Nn * ...
            (-YthetaExpf*((1i)^(nn+1)*bidx + (-1i)^(nn+1)*qidx) ...
            + YphiExpf*((1i)^nn*aidx + (-1i)^nn*pidx) );
        end
      end

      E=[zeros(size(Etheta)),Etheta,Ephi].';
      H=[zeros(size(Htheta)),Htheta,Hphi].';
      rtp = [zeros(size(Etheta)), theta(:), phi(:)].';

      % SI-ify units of H
      H = H * -1i;

      % Package output
      E = ott.utils.FieldVector(rtp, E, 'spherical');
      H = ott.utils.FieldVector(rtp, H, 'spherical');
    end

    function [E, H, data] = ehfieldRtp(beam, rtp, varargin)
      % Calculates the E and H field at specified locations
      %
      % Usage
      %   [E, H] = beam.ehfieldRtp(rtp, ...) calculates the complex field
      %   at locations xyz (3xN matrix of spherical coordinates).
      %   Returns 3xN matrices for the E and H field at these locations.
      %
      %   [E, H, data] = beam.emFieldXyz(...) returns data used in the
      %   calculation.  See 'data' optional parameter.
      %
      % Optional named arguments:
      %   - calcE   bool   calculate E field (default: true)
      %
      %   - calcH   bool   calculate H field (default: nargout == 2)
      %
      %   - saveData bool  save data for repeated calculation (default: false)
      %
      %   - data    data   data saved for repeated calculation.
      %
      %   - coord   str    coordinates to use for calculated field
      %    'cartesian' (default) or 'spherical'
      %
      % If either calcH or calcE is false, the function still returns
      % E and H as matrices of all zeros for the corresponding field.
      %
      % If internal fields are calculated only the theta and phi components
      % of E are continuous at the boundary. Conversely, only the kr
      % component of D is continuous at the boundary.
      %
      % This function is based on the emField function from OTTv1.

      p = inputParser;
      p.addParameter('calcE', true);
      p.addParameter('calcH', nargout >= 2);
      p.addParameter('saveData', false);
      p.addParameter('data', []);
      p.parse(varargin{:});

      % Scale the locations by the wave number (unitless coordinates)
      rtp(1, :) = rtp(1, :) * abs(beam.wavenumber);

      % Get the indices required for the calculation
      [n,m]=ott.utils.combined_index(find(any(abs(beam.a)|abs(beam.b), 2)));

      ci = ott.utils.combined_index(n, m);
      [a, b] = beam.getCoefficients(ci);

      % Coherently combine coefficients before calculation
      if strcmpi(beam.array_type, 'coherent')
        a = sum(a, 2);
        b = sum(b, 2);
      end

      % Get unique rtp for faster calculation
      [r_new,~,indR]=unique(rtp(1, :).');
      [theta_new,~,indTheta]=unique(rtp(2, :).');
      [phi_new,~,indPhi]=unique(rtp(3, :).');

      % Remove zeros
      r_new(r_new == 0) = 1e-15;

      % Allocate memory for outputs
      E = {zeros(size(rtp)).'};
      H = {zeros(size(rtp)).'};

      if numel(beam) > 1
        E = repmat(E, size(beam));
        H = repmat(H, size(beam));
      end

      un = unique(n);

      % Allocate memory for output data
      data = [];
      if p.Results.saveData
        data = zeros(numel(indTheta), 0);
      end

      % Start a counter for accessing the data
      if ~isempty(p.Results.data)
        dataCount = 0;
      end

      for nn = 1:max(un)
        Nn = 1/sqrt(nn*(nn+1));
        vv=find(n==nn);

        if ~isempty(vv)

          kr=r_new(indR);

          if isempty(p.Results.data)

            [Y,Ytheta,Yphi] = ott.utils.spharm(nn,m(vv),theta_new, ...
                zeros(size(theta_new)));

            switch beam.basis
              case 'incoming'
                [hn,dhn]=ott.utils.sbesselh2(nn,r_new);
                hn = hn ./ 2;
                dhn = dhn ./ 2;

              case 'outgoing'
                [hn,dhn]=ott.utils.sbesselh1(nn,r_new);
                hn = hn ./ 2;
                dhn = dhn ./ 2;

              case 'regular'
                [hn,dhn]=ott.utils.sbesselj(nn,r_new);

              otherwise
                error('Unknown beam type');
            end

            [M,PHI]=meshgrid(1i*m(vv),phi_new);

            expimphi=exp(M.*PHI);

            hnU=hn(indR);
            dhnU=dhn(indR);

            [jnr,djnr]=ott.utils.sbesselj(nn,r_new);
            jnrU=jnr(indR);
            djnrU=djnr(indR);

            % Create full Y, Ytheta, Yphi, expimphi matrices (opt, R2018a)
            expimphif = expimphi(indPhi, :);
            YExpf = Y(indTheta, :).*expimphif;
            YthetaExpf = Ytheta(indTheta, :).*expimphif;
            YphiExpf = Yphi(indTheta, :).*expimphif;

            % Save the data if requested
            if p.Results.saveData
              data(:, end+1) = hnU;
              data(:, end+1) = dhnU;
              data(:, end+(1:size(Ytheta, 2))) = YExpf;
              data(:, end+(1:size(Ytheta, 2))) = YthetaExpf;
              data(:, end+(1:size(Ytheta, 2))) = YphiExpf;
            end

          else

            % Load the data if present
            hnU = p.Results.data(:, dataCount+1);
            dataCount = dataCount + 1;
            dhnU = p.Results.data(:, dataCount+1);
            dataCount = dataCount + 1;
            YExpf = p.Results.data(:, dataCount+(1:length(vv)));
            dataCount = dataCount + length(vv);
            YthetaExpf = p.Results.data(:, dataCount+(1:length(vv)));
            dataCount = dataCount + length(vv);
            YphiExpf = p.Results.data(:, dataCount+(1:length(vv)));
            dataCount = dataCount + length(vv);

          end

          for ii = 1:size(a, 2)
            pidx = full(a(vv, ii));
            qidx = full(b(vv, ii));

            % Now we use full matrices, we can use matmul (opt, R2018a)
            if p.Results.calcE
              E{ii}(:,1)=E{ii}(:,1)+Nn*nn*(nn+1)./kr.*hnU.*YExpf*qidx(:);
              E{ii}(:,2)=E{ii}(:,2)+Nn*(hnU.*YphiExpf*pidx(:) ...
                  + dhnU.*YthetaExpf*qidx(:));
              E{ii}(:,3)=E{ii}(:,3)+Nn*(-hnU.*YthetaExpf*pidx(:) ...
                  + dhnU.*YphiExpf*qidx(:));
            end

            if p.Results.calcH
              H{ii}(:,1)=H{ii}(:,1)+Nn*nn*(nn+1)./kr.*hnU.*YExpf*pidx(:);
              H{ii}(:,2)=H{ii}(:,2)+Nn*((hnU(:).*YphiExpf)*qidx(:) ...
                  +(dhnU(:).*YthetaExpf)*pidx(:));
              H{ii}(:,3)=H{ii}(:,3)+Nn*((-hnU(:).*YthetaExpf)*qidx(:) ...
                  +(dhnU(:).*YphiExpf)*pidx(:));
            end
          end
        end
      end

      for ii = 1:size(a, 2)
        H{ii} = H{ii}.';
        E{ii} = E{ii}.';

        H{ii}=-1i*H{ii}; %LOOK HERE TO FIX STUFF

        % Package output
        E{ii} = ott.utils.FieldVector(rtp, E{ii}, 'spherical');
        H{ii} = ott.utils.FieldVector(rtp, H{ii}, 'spherical');
      end

      % Discard cell array if we don't need it
      if size(a, 2) == 1
        E = E{1};
        H = H{1};
      end
    end

    function [E, H, data] = ehfield(beam, xyz, varargin)
      % Calculates the E and H field at specified locations
      %
      % Usage
      %   [E, H] = beam.emFieldXyz(xyz, ...) calculates the complex field
      %   at locations xyz (3xN matrix of Cartesian coordinates).
      %   Returns 3xN :class:`ott.utils.FieldVector` for the E and H fields.
      %
      %   [E, H, data] = beam.emFieldXyz(...) returns data used in the
      %   calculation.  See 'data' optional parameter.
      %
      % For further details, see :meth:`ehfieldRtp`.

      rtp = ott.utils.xyz2rtp(xyz);
      [E, H, data] = beam.ehfieldRtp(rtp, varargin{:});
    end

    function varargout = visualise(beam, varargin)
      % Create a visualisation of the beam
      %
      % Usage
      %   beam.visualise(...) displays an image of the beam in the current
      %   figure window.
      %
      %   im = beam.visualise(...) returns a image of the beam.
      %   If the beam object is an array, returns an image for each beam.
      %
      % Optional named arguments
      %   - size (2 numeric) -- Number of rows and columns in image.
      %     Default: ``[80, 80]``.
      %
      %   - field (enum) -- Type of visualisation to generate, see
      %     :meth:`VisualisationData` for valid options.
      %     Default: ``'irradiance'``.
      %
      %   - axis (enum|cell) -- Describes the slice to visualise.
      %     Can either be 'x', 'y' or 'z' for a plane perpendicular to
      %     the x, y and z axes respectively; or a cell array with 2
      %     or 3 unit vectors (3x1) for x, y, [z] directions.
      %     Default: ``'z'``.
      %
      %   - offset (numeric) -- Plane offset along axis (default: 0.0)
      %
      %   - range (numeric|cell) -- Range of points to visualise.
      %     Can either be a cell array { x, y }, two scalars for
      %     range [-x, x], [-y, y] or 4 scalars [ x0, x1, y0, y1 ].
      %     Default: ``[1, 1].*nmax2ka(Nmax)/abs(wavenumber)``.
      %
      %   - mask (function_handle) Describes regions to remove from the
      %     generated image.  Function should take one argument for the
      %     3xN field xyz field locations and return a logical array mask.
      %     Default: ``[]``.
      %
      %   - axes (axes handle) -- Axes to place the visualisation in.
      %     If empty, no visualisation is generated.
      %     Default: ``gca()`` if ``nargout == 0`` otherwise ``[]``.
      %
      %   - combine (enum|empty) -- Method to use when combining beams.
      %     Can either be emtpy (default), 'coherent' or 'incoherent'.

      import ott.utils.nmax2ka;

      % Default parameter changes for this function
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('range', [1,1].*nmax2ka(max([beam.Nmax])) ...
          ./abs(min([beam.wavenumber])));
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      [varargout{1:nargout}] = visualise@ott.beam.Beam(...
          beam, unmatched{:}, 'range', p.Results.range);
    end

    function varargout = visualiseFarfield(beam, varargin)
      % Create a visualisation of the beam by projecting the far-field
      % onto a plane.
      %
      % Usage
      %   beam.visualiseFarfield(...) displays an image of the beam
      %   in the current axes.
      %
      % Optional named arguments
      %   - basis (enum) -- Basis to use for visualisation.  Can either
      %     be 'incoming' or 'outgoing'.  Defaults to the beam basis
      %     unless 'regular', in which case uses 'incoming' basis.
      %
      % See :class:`ott.beam.Beam` for usage/parameters.

      [beam, unmatched] = beam.farfieldVisHelper(varargin{:});
      [varargout{1:nargout}] = visualiseFarfield@ott.beam.Beam(...
          beam, unmatched{:});
    end

    function varargout = visualiseFarfieldSphere(beam, varargin)
      % Generate a spherical surface visualisation of the far-field
      %
      % Usage
      %   beam.visualiseFarfieldSphere(...)
      %   Generate a visualisation of the far-field in the current axes.
      %
      % Optional named arguments
      %   - basis (enum) -- Basis to use for visualisation.  Can either
      %     be 'incoming' or 'outgoing'.  Defaults to the beam basis
      %     unless 'regular', in which case uses 'incoming' basis.
      %
      % See :class:`ott.beam.Beam` for usage/parameters.

      [beam, unmatched] = beam.farfieldVisHelper(varargin{:});
      [varargout{1:nargout}] = visualiseFarfieldSphere@ott.beam.Beam(...
          beam, unmatched{:});
    end

    function varargout = visualiseFarfieldSlice(beam, phi, varargin)
      % Generate a 2-D slice through the far-field
      %
      % Usage
      %   beam.visualiseFarfieldSlice(phi, ...)
      %   Generates a 2-D slice at angle phi around the z-axis.
      %   Plots into the current axes.
      %
      % Optional named arguments
      %   - basis (enum) -- Basis to use for visualisation.  Can either
      %     be 'incoming' or 'outgoing'.  Defaults to the beam basis
      %     unless 'regular', in which case uses 'incoming' basis.
      %
      % See :class:`ott.beam.Beam` for usage/parameters.

      [beam, unmatched] = beam.farfieldVisHelper(varargin{:});
      [varargout{1:nargout}] = visualiseFarfieldSlice@ott.beam.Beam(...
          beam, unmatched{:});
    end

    function [moments, ints] = intensityMoment(beam, varargin)
      % Calculate moment of the beam intensity in the far-field
      %
      % Usage
      %   [moment, int] = beam.intensityMoment(...)
      %
      % Optional named arguments
      %   - basis (enum) -- Basis to use for visualisation.  Can either
      %     be 'incoming' or 'outgoing'.  Defaults to the beam basis
      %     unless 'regular', in which case uses 'incoming' basis.
      %
      % See :class:`ott.beam.Beam` for parameters/usage.

      [beam, unmatched] = beam.farfieldVisHelper(varargin{:});
      [varargout{1:nargout}] = intensityMoment@ott.beam.Beam(...
          beam, unmatched{:});
    end

    % TODO: Ensure sparsity is respected
    % TODO: Add a getDenseBeamVectors function (add to getCoefficients?)

    function b = issparse(beam)
      % Returns true if the VSWF data is sparse
      %
      % Usage
      %   b = issparse(beam)

      b = issparse(beam.getCoefficients());
    end

    function beam = full(beam)
      % Convert the VSWF data to a full matrix
      %
      % Usage
      %   beam = full(beam)

      beam = beam.setCoefficients(full(beam.getCoefficients()));
    end

    function beam = sparse(beam)
      % Convert the VSWF data to a sparse matrix
      %
      % This function doesn't change the data.  For a method that removes
      % near-zeros elements, see :meth:`makeSparse`.
      %
      % Usage
      %   beam = sparse(beam)

      beam = beam.setCoefficients(sparse(beam.getCoefficients()));
    end

    function beam = makeSparse(beam, varargin)
      % Make the beam sparse by removing near-zero power elements
      %
      % Usage
      %   beam = beam.makeSparse(...)
      %
      % Optional named arguments
      %   - tolerance (numeric) -- Tolerance for pruning elements.
      %     Default: ``1.0e-15``.

      ott.utils.nargoutCheck(beam, nargout);

      p = inputParser;
      p.addParameter('tolerance', 1.0e-15, @isnumeric);
      p.parse(varargin{:});

      pw = abs(beam.a).^2 + abs(beam.b).^2;
      non_zero = pw > p.Results.tolerance*max(pw);

      beam.a(~non_zero) = 0;
      beam.b(~non_zero) = 0;

      beam.a = sparse(beam.a);
      beam.b = sparse(beam.b);
    end

    function nbeam = shrinkNmax(beam, varargin)
      % Reduces the size of the beam while preserving power
      %
      % Usage
      %   beam = beam.shrinkNmax(...)
      %
      % Optional named arguments
      %   - tolerance (numeric) -- Tolerance to use for power loss.
      %     Default: ``1.0e-6``.

      ott.utils.nargoutCheck(beam, nargout);

      p = inputParser;
      p.addParameter('tolerance', 1.0e-6, @isnumeric);
      p.parse(varargin{:});

      amagA = full(sum(sum(abs(beam.a).^2)));
      bmagA = full(sum(sum(abs(beam.b).^2)));

      for ii = 1:beam.Nmax

        total_orders = ott.utils.combined_index(ii, ii);
        nbeam = beam;
        nbeam.a = nbeam.a(1:total_orders);
        nbeam.b = nbeam.b(1:total_orders);

        amagB = full(sum(sum(abs(nbeam.a).^2)));
        bmagB = full(sum(sum(abs(nbeam.b).^2)));

        aapparent_error = abs( amagA - amagB )/amagA;
        bapparent_error = abs( bmagA - bmagB )/bmagA;

        if aapparent_error < p.Results.tolerance && ...
            bapparent_error < p.Results.tolerance
          break;
        end
      end
    end

    function beam = setNmax(beam, nmax, varargin)
      % Resize the beam, with additional options
      %
      % Usage
      %   beam = beam.setNmax(nmax, ...)   or    beam.Nmax = nmax
      %   Set the Nmax, a warning is issued if truncation occurs.
      %
      % Optional named arguments
      %   - tolerance (numeric) -- Specify the tolerance to use for
      %     power loss warnings.  Default: ``1.0e-6``.
      %
      %   - powerloss (enum) -- Action to take when beam power is lost.
      %     Can be one of 'ignore', 'warn' or 'error'.
      %     Default: ``'warn'``.

      p = inputParser;
      p.addParameter('tolerance', 1.0e-6);
      p.addParameter('powerloss', 'warn');
      p.parse(varargin{:});

      total_orders = ott.utils.combined_index(nmax, nmax);
      if size(beam.a, 1) > total_orders

        amagA = full(sum(sum(abs(beam.a).^2)));
        bmagA = full(sum(sum(abs(beam.b).^2)));

        beam.a = beam.a(1:total_orders, :);
        beam.b = beam.b(1:total_orders, :);

        amagB = full(sum(sum(abs(beam.a).^2)));
        bmagB = full(sum(sum(abs(beam.b).^2)));

        if ~strcmpi(p.Results.powerloss, 'ignore')

          aapparent_error = abs( amagA - amagB )/amagA;
          bapparent_error = abs( bmagA - bmagB )/bmagA;

          if aapparent_error > p.Results.tolerance || ...
              bapparent_error > p.Results.tolerance
            if strcmpi(p.Results.powerloss, 'warn')
              warning('ott:beam:vswf:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            elseif strcmpi(p.Results.powerloss, 'error')
              error('ott:beam:vswf:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            else
              error('ott:beam:vswf:Bsc:setNmax:truncation', ...
                'powerloss should be one of ignore, warn or error');
            end
          end
        end
      elseif size(beam.a, 1) < total_orders
        [arow_index,acol_index,aa] = find(beam.a);
        [brow_index,bcol_index,ba] = find(beam.b);
        beam.a = sparse(arow_index,acol_index,aa,total_orders,numel(beam));
        beam.b = sparse(brow_index,bcol_index,ba,total_orders,numel(beam));
      end
    end

    function bsc = applyTransformation(bsc, varargin)
      % Apply a transformation to the beam.
      %
      % Rotation is applied first, followed by translation.
      %
      % Units for the coordinates should be consistent with the
      % beam wave number (i.e., if the beam was created by specifying
      % wavelength in units of meters, distances here should also be
      % in units of meters).
      %
      % By default, this operation demotes the beam type to vswf.Bsc if a
      % transformation is applied.  Otherwise the beam is unchanged.
      % The position and rotation are set to itendity after transformation.
      %
      % Usage
      %   bsc = bsc.applyTransformation(...)
      %
      % Optional named arguments
      %   - rotation (3x3N numeric) -- Rotation to apply.
      %     Default: ``bsc.position``.
      %
      %   - position (3xN numeric) -- Translation to apply.
      %     Default: ``bsc.rotation``.

      p = inputParser;
      p.addParameter('position', bsc.position);
      p.addParameter('rotation', bsc.rotation);
      p.parse(varargin{:});

      ott.utils.nargoutCheck(bsc, nargout);

      % TODO: Pre-rotation

      bsc = bsc.applyRotation(p.Results.rotation);
      bsc = bsc.applyTranslation(p.Results.position);
    end

    function [bsc, D] = applyRotation(bsc, varargin)
      % Apply rotation to the beam shape coefficients.
      %
      % This operation demotes the beam type to vswf.Bsc if a
      % rotation is applied.  Otherwise the beam is unchanged.
      %
      % Usage
      %   bsc = bsc.applyRotation(...)
      %   Applies ``bsc.rotation`` to the beam shape coefficients.
      %
      %   bsc = bsc.applyRotation(R, ...)
      %   Applies a specific rotation.  R should be a 3x3N matrix.
      %
      %   [bsc, D] = bsc.applyRotation(...)
      %   Additionally returns the wigner rotation matrix used for
      %   the operation.  Note, is generally faster *not* to use the
      %   wigner rotation matrix for inverse rotations.
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Requested minimum Nmax for rotated beam.
      %     Default: ``bsc.Nmax``.
      %
      %   - wigner (numeric|cell) -- A wigner rotation matrix or
      %     a cell array of wigner rotation matrices to use for the
      %     rotation.

      p = inputParser;
      p.addOptional('R', [], @isnumeric);
      p.addParameter('wigner', []);
      p.addParameter('Nmax', bsc.Nmax);
      p.parse(varargin{:});

      ott.utils.nargoutCheck(bsc, nargout);

      if isempty(p.Results.wigner)

        R = p.Results.R;
        if isempty(R)
          R = bsc.rotation;
        else
          assert(isnumeric(R) && ismatrix(R) && size(R, 1) == 3 ...
              && mod(size(R, 2), 3) == 0, ...
              'R must be 3x3N matrix');
        end

        Nrots = size(R, 2)/3;
        Nbeams = size(bsc, 2);

        assert(Nrots == 1 || Nbeams == 1 || Nrots == Nbeams, ...
            'Number of rotations must match number of beams or be scalar');

        Nwork = max([Nrots, Nbeams]);
        if Nwork > 1
          D = cell(1, Nwork);
          for ii = 1:Nwork
            D = ott.utils.wigner_rotation_matrix(...
                max([bsc.Nmax, p.Results.Nmax]), R(:, (1:3) + (ii-1)*3));
          end
        else
          D = ott.utils.wigner_rotation_matrix(...
              max([bsc.Nmax, p.Results.Nmax]), R);
        end
      elseif ~isempty(p.Results.R)
        error('Too many arguments: Must only supply wigner or R, not both');
      else
        D = p.Results.wigner;

        if iscell(D)
          Nrots = numel(D);
        else
          Nrots = 1;
        end
        Nbeams = size(bsc, 2);
        Nwork = max([Nrots, Nbeams]);

        assert(Nrots == 1 || Nbeams == 1 || Nrots == Nbeams, ...
            'Number of rotations must match number of beams or be scalar');
      end

      % Apply wigner rotation matrices
      if iscell(D)

        % Check for work to do
        hasWork = false;
        for ii = 1:numel(D)
          hasWork = hasWork | any(any(D ~= eye(size(D))));
        end

        if hasWork
          ibsc = bsc;
          bsc = ott.beam.vswf.Bsc.empty([1, numel(D)]);
          sz = size(bsc.a, 1);

          for ii = 1:numel(D)
            bsc(ii) = D{ii}(:, 1:sz) * ibsc;
          end
        else
          % Nothing to do (leave type unchanged)
        end
      else
        if ~all(all(D == eye(size(D))))
          sz = size(bsc.a, 1);
          bsc = ott.beam.vswf.Bsc(D(:, 1:sz) * bsc);
        else
          % Nothing to do (leave type unchanged)
        end
      end

      % Set the rotation to identity
      bsc.rotation = eye(3);
    end

    function [bsc, Az, Bz, D] = applyTranslation(bsc, varargin)
      % Apply translation to the beam shape coefficients.
      %
      % Applies rotations by first rotating the beam z-axis to align
      % to the translation axis, then translating along z, before
      % rotating back to the original orientation.
      %
      % Units for the coordinates should be consistent with the
      % beam wave number (i.e., if the beam was created by specifying
      % wavelength in units of meters, distances here should also be
      % in units of meters).
      %
      % This operation demotes the beam type to vswf.Bsc if a
      % translation is applied.  Otherwise the beam is unchanged.
      %
      % Usage
      %   bsc = bsc.applyTranslation(...)
      %   Applies ``bsc.position`` to the beam shape coefficients.
      %
      %   bsc = bsc.applyTranslation(P, ...)
      %   Applies a specific rotation.  P should be a 3xN matrix.
      %
      %   [bsc, Az, Bz, D] = bsc.applyRotation(...)
      %   Additionally returns the translation matrices.
      %
      % Optional named arguments
      %   - Nmax (numeric) -- Requested minimum Nmax for translated beam.
      %     The Nmax limit is applied during the translation. The first
      %     rotation step uses the full Nmax, the last uses the new Nmax.
      %     Ignored when multiple outputs are requested.
      %     Default: ``bsc.Nmax``.
      %
      %   - AB (2xN|3xN cell) -- Cell array of translation matrices to
      %     apply.  For 2xN cell array, assumes {A; B}, for 3xN cell
      %     array {Az; Bz; D} where D is the rotation to the z axis.

      p = inputParser;
      p.addOptional('P', []);
      p.addParameter('Nmax', bsc.Nmax);
      p.addParameter('AB', {});
      p.parse(varargin{:});

      ott.utils.nargoutCheck(bsc, nargout);

      if isempty(p.Results.AB)

        P = p.Results.P;
        if isempty(P)
          P = bsc.position;
        else
          assert(isnumeric(P) && ismatrix(P) && size(P, 1) == 3, ...
              'P must be 3xN numeric matrix');
        end

        Nrots = size(P, 2);
        Nbeams = size(bsc, 2);
        Nwork = max([Nrots, Nbeams]);

        assert(Nrots == 1 || Nbeams == 1 || Nrots == Nbeams, ...
            'Number of positions must match number of beams or be scalar');

        Nmax = p.Results.Nmax;
        if nargout > 1
          Nmax = bsc.Nmax;
        end

        % TODO: Optimisation: Unique rotation directions

        % Check for work to do
        hasWork = false;
        for ii = 1:numel(Nrots)
          hasWork = hasWork | any(P(:, ii) ~= [0;0;0]);
        end

        Az = cell(1, Nrots);
        Bz = cell(1, Nrots);
        D = cell(1, Nrots);

        if hasWork

          ibsc = bsc;
          bsc = ott.beam.vswf.Bsc.empty([1, Nrots]);

          for ii = 1:Nrots

            rtp = ott.utils.xyz2rtp(P(:, ii));
            R = ott.utils.rotz(rtp(3)*180/pi) * ott.utils.roty(rtp(2)*180/pi);

            if nargout > 1
              [newbeam, D{ii}] = ibsc.applyRotation(R);
              [newbeam, Az{ii}, Bz{ii}] = newbeam.applyZTranslation(...
                  rtp(1), 'Nmax', Nmax);
              newbeam = D{ii}' * newbeam;
            else
              newbeam = ibsc.applyRotation(R);
              newbeam = newbeam.applyZTranslation(rtp(1), 'Nmax', Nmax);
              newbeam = newbeam.applyRotation(R.');
            end

            bsc(ii) = newbeam;
          end
        else
          bsc = repmat(bsc, 1, Nrots);
          Az = repmat({1}, 1, Nrots);
          Bz = repmat({0}, 1, Nrots);
          D = repmat({1}, 1, Nrots);
        end

        if Nrots == 1 && nargout > 1
          Az = Az{1};
          Bz = Bz{1};
          D = D{1};
        end

      elseif ~isempty(p.Results.P)
        error('Too many arguments: only one of P or AB must be supplied');
      else

        AB = p.Results.AB;
        assert(iscell(AB) && (size(AB, 1) == 2 || size(AB, 1) == 3), ...
          'AB must be 2xN or 3xN cell array');

        % Unpack AB coefficients
        A = AB(1, :);
        B = AB(2, :);

        Nrots = numel(A);
        Nbeams = size(bsc, 2);
        Nwork = max([Nrots, Nbeams]);

        assert(Nrots == 1 || Nbeams == 1 || Nrots == Nbeams, ...
            'Number of AB must match number of beams or be scalar');

        ibsc = bsc;
        bsc = ott.beam.vswf.Bsc.empty([1, numel(A)]);
        sz = size(ibsc.a, 1);

        if size(AB, 1) == 2
          for ii = 1:numel(A)

            A2 = A{ii}(:, 1:sz);
            B2 = B{ii}(:, 1:sz);

            bsc(ii) = [A2 B2; B2 A2] * ibsc;
          end
        else
          D = AB(3, :);
          for ii = 1:numel(A)

            % Need all same size for rotation
            A2 = A{ii}(1:sz, 1:sz);
            B2 = B{ii}(1:sz, 1:sz);
            D2 = D{ii}(1:sz, 1:sz);

            bsc = [bsc, D2' * ([A2 B2; B2 A2] * (D2 * ibsc))];
          end
        end
      end

      % Set the position to 0
      bsc.position = [0;0;0];
    end

    function [bsc, A, B] = applyZTranslation(bsc, z, varargin)
      % Apply a translation along the beam z-axis.
      %
      % Units for the coordinates should be consistent with the
      % beam wave number (i.e., if the beam was created by specifying
      % wavelength in units of meters, distances here should also be
      % in units of meters).
      %
      % Usage
      %   [bsc, A, B] = bsc.applyZTranslation(z, ...)
      %
      %   [bsc, AB] = bsc.applyZTranslation(z, ...)
      %   Returns AB matrices packed for convenient application.
      %
      % Optional named parameters
      %   - Nmax (numeric) -- Requested minimum Nmax for translated beam.
      %     Default: ``bsc.Nmax``.

      p = inputParser;
      p.addParameter('Nmax', bsc.Nmax);
      p.parse(varargin{:});

      ott.utils.nargoutCheck(bsc, nargout);

      % Check amount of work
      Nz = numel(z);
      Nbeams = numel(bsc);
      assert(Nz == 1 || Nbeams == 1 || Nz == Nbeams, ...
          'Number of positions and beams must be equal or scalar');

      ibsc = bsc;
      bsc = ott.beam.vswf.Bsc.empty([1, numel(z)]);

      % Add a warning when the beam is translated outside nmax2ka(Nmax)
      % The first time may be OK, the second time does not have enough
      % information.
      if any(ibsc.absdz > ott.utils.nmax2ka(ibsc.Nmax)./ibsc.wavenumber)
        warning('ott:beam:vswf:Bsc:translateZ:outside_nmax', ...
            'Repeated translation of beam outside Nmax region');
      end
      bsc.absdz = ibsc.absdz + abs(z);

      % Convert to beam units
      z = z * ibsc.wavenumber / 2 / pi;

      A = cell(1, Nz);
      B = cell(1, Nz);

      for ii = 1:numel(z)

        [A{ii}, B{ii}] = ibsc.applyZTranslationInternal(...
            z(ii), [p.Results.Nmax, ibsc.Nmax]);
        bsc(ii) = [A{ii}, B{ii}; B{ii}, A{ii}] * ibsc;

      end
      
      % Pack output for convinient application
      if nargout < 3
        A = cellfun(@(a, b) [a b; b a], A, B, 'UniformOutput', false);
      end
      
      if numel(z) == 1
        A = A{1};
        B = B{1};
      end
    end

    function beam = setCoefficients(beam, a, b)
      % Set the `a` and `b` coefficients
      %
      % Usage
      %   beam = beam.setCoefficients(a, b)
      %
      %   beam = beam.setCoefficients(ab)
      %   Set coefficients from a [a; b] vector.
      %
      %   beam = beam.setCoefficients(other_beam)
      %   Get the beam coefficients from another beam
      
      ott.utils.nargoutCheck(beam, nargout);
      
      if nargin == 2
        if isa(a, 'ott.beam.vswf.Bsc')
          b = a.b;
          a = a.a;
        else
          assert(isnumeric(a) && ismatrix(a), 'ab must be numeric matrix');
          assert(mod(size(a, 1), 2) == 0, 'ab must be 2NxM in size');
          b = a(end/2+1:end, :);
          a = a(1:end/2, :);
        end
      end

      assert(all(size(a) == size(b)), 'size of a and b must match');
      assert(nargout == 1, 'Expected one output from function');

      assert(isempty(a) || ...
          (size(a, 1) >= 3 && ...
          sqrt(size(a, 1)+1) == floor(sqrt(size(a, 1)+1))), ...
        'number of multipole terms must be empty, 3, 8, 15, 24, ...');

      beam.a = a;
      beam.b = b;
    end

    function varargout = getCoefficients(beam, ci)
      % Gets the beam a/b coefficients
      %
      % Usage
      %   ab = beam.getCoefficients() gets the beam coefficients packed
      %   into a single vector, suitable for multiplying by a T-matrix.
      %
      %   [a, b] = beam.getCoefficients() get the coefficients in two
      %   beam vectors.
      %
      %   [...] = beam.getCoefficients(ci) behaves as above but only returns
      %   the requested beam cofficients a(ci) and b(ci).

      % If ci omitted, return all a and b
      if nargin == 1
        ci = 1:size(beam.a, 1);
      end

      a = beam.a(ci, :);
      b = beam.b(ci, :);

      if nargout == 1
        varargout{1} = [a; b];
      else
        [varargout{1:2}] = deal(a, b);
      end
    end

    function [n, m] = getModeIndices(beam)
      %GETMODEINDICES gets the mode indices
      [n, m] = ott.utils.combined_index([1:size(beam.a, 1)].');
      if nargout == 1
        n = [n; m];
      end
    end

    function beam = mrdivide(beam,o)
      %MRDIVIDE (op) divide the beam coefficients by a scalar
      beam.a = beam.a / o;
      beam.b = beam.b / o;
    end

    function beam = mtimes(a,b)
      % Provides beam matrix and scalar multiplication
      %
      % Usage
      %   beam = scalar * beam   or   beam = beam * scalar
      %   Scalar multiplication of beam shape coefficients.
      %
      %   beam = matrix * beam
      %   Matrix multiplication.  Matrix can either have the same
      %   number of columns as the `a` and `b` beam shape coefficients
      %   or half as many rows.  In other words, the resulting beam
      %   is either `[M*a; M*b]` or `M*[a;b]`.

      if isa(a, 'ott.beam.vswf.Bsc')
        beam = a;
        beam.a = beam.a * b;
        beam.b = beam.b * b;
      else
        beam = b;
        if size(a, 2) == 2*size(beam.a, 1)
          ab = a * [beam.a; beam.b];
          beam.a = ab(1:size(ab, 1)/2, :);
          beam.b = ab(1+size(ab, 1)/2:end, :);
        else
          beam.a = a * beam.a;
          beam.b = a * beam.b;
        end
      end
    end

    function beam = times(beam, rv)
      % Provides element-wise multiplication of BSC coefficients.
      %
      % Usage
      %   beam = beam .* row_vec
      %   Multiplies each beam by the elements of `row_vec`.

      assert(isa(beam, 'ott.beam.vswf.Bsc'), 'first argument must be Bsc');
      assert(isnumeric(rv) && ismatrix(rv) && size(rv, 1) == 1 ...
          && size(rv, 2) == numel(beam), ...
          'second argument must be row vector matching N-beams');

      beam.a = beam.a .* rv;
      beam.b = beam.b .* rv;
    end

    function beam = minus(beam1, beam2)
      %MINUS subtract two beams

      if beam1.Nmax > beam2.Nmax
        beam2.Nmax = beam1.Nmax;
      elseif beam2.Nmax > beam1.Nmax
        beam1.Nmax = beam2.Nmax;
      end

      beam = beam1;
      beam.a = beam.a - beam2.a;
      beam.b = beam.b - beam2.b;
    end

    function beam = sum(beamin, dim)
      % Sum beam coefficients
      %
      % Usage
      %   beam = sum(beam)
      %
      %   beam = beam.sum()
      %
      %   beam = sum([beam1, beam2, ...], dim) sums the given beams,
      %   similar to Matlab's ``sum`` builtin.  ``dim`` is the dimension
      %   to sum over (optional).

      if numel(beamin) > 1
        % beam is an array

        % Handle default value for dimension
        if nargin < 2
          if isvector(beamin)
            dim = find(size(beamin) > 1, 1);
          elseif ismatrix(beamin) == 2
            dim = 2;
          else
            dim = find(size(beamin) > 1, 1);
          end
        end

        % Select the first row in our dimension
        subs = [repmat({':'},1,dim-1), 1, ...
          repmat({':'},1,ndims(beamin)-dim)];
        S = struct('type', '()', 'subs', {subs});
        beam = subsref(beamin, S);

        % Add each beam
        for ii = 2:size(beamin, dim)
        subs = [repmat({':'},1,dim-1), ii, ...
          repmat({':'},1,ndims(beamin)-dim)];
          S = struct('type', '()', 'subs', {subs});
          beam = beam + subsref(beamin, S);
        end

      else
        % Beam union
        beam = beamin;
        beam.a = sum(beam.a, 2);
        beam.b = sum(beam.b, 2);
      end
    end

    function beam = plus(b1, b2)
      % Add two beams together
      %
      % If the beams are both Bsc beams, adds the field coefficients.
      % Otherwise, deffers to the base class plus operation.

      if isa(b1, 'ott.beam.vswf.Bsc') && isa(b2, 'ott.beam.vswf.Bsc')
        beam = plusInternal(b1, b2);
      else
        beam = plus@ott.beam.utils.ArrayType(b1, b2);
      end
    end
  end

  methods (Hidden)
    function bsc = translateXyzInternal(bsc, xyz, varargin)
      % Apply translations to position property.
      %
      % If we have more than one beam or translation, converts the
      % Bsc object to an array and applies translations to array elements.

      if size(xyz, 2) > 1
        bsc = ott.beam.vswf.Array(bsc);
        bsc = bsc.translateXyz(xyz, varargin{:});
      else
        bsc = translateXyzInternal@ott.beam.utils.ArrayType(...
            bsc, xyz, varargin{:});
      end
    end

    function [A, B] = applyZTranslationInternal(beam, z, Nmax)
      % Generate the translation matrices used by translate z.
      %
      % This behaviour can be different with some beams.
      %
      % Units for the coordinates should be consistent with the
      % beam wave number (i.e., if the beam was created by specifying
      % wavelength in units of meters, distances here should also be
      % in units of meters).
      %
      % Usage
      %   [A, B] = translateZ_type_helper(beam, z, Nmax) calculates the
      %   translation matrices for distance z with Nmax
      %
      % Usage may change in future releases.

      % Determine beam type
      switch beam.basis
        case 'incoming'
          translation_type = 'sbesselh2';
        case 'outgoing'
          translation_type = 'sbesselh1';
        case 'regular'
          translation_type = 'sbesselj';
      end

      % Calculate tranlsation matrices
      [A, B] = ott.utils.translate_z(Nmax, z, 'type', translation_type);

    end

    function beam = catInternal(dim, beam, varargin)
      % Concatenate beams

      assert(dim == 2, 'Only horzcat (dim=2) supported for now');

      other_a = {};
      other_b = {};
      for ii = 1:length(varargin)
        other_a{ii} = varargin{ii}.a;
        other_b{ii} = varargin{ii}.b;
      end

      beam.a = cat(dim, beam.a, other_a{:});
      beam.b = cat(dim, beam.b, other_b{:});
    end

    function beam = plusInternal(beam1, beam2)
      %PLUS add two beams together

      if beam1.Nmax > beam2.Nmax
        beam2.Nmax = beam1.Nmax;
      elseif beam2.Nmax > beam1.Nmax
        beam1.Nmax = beam2.Nmax;
      end

      beam = beam1;
      beam.a = beam.a + beam2.a;
      beam.b = beam.b + beam2.b;
    end

    function beam = subsrefInternal(beam, subs)
      % Get the subscripted beam

      if numel(subs) > 1
        if subs{1} == 1 || strcmpi(subs{1}, ':')
          subs = subs(2:end);
        end
        assert(numel(subs) == 1, 'Only 1-D indexing supported for now');
      end

      beam.a = beam.a(:, subs{:});
      beam.b = beam.b(:, subs{:});
    end

    function beam = subsasgnInternal(beam, subs, rem, other)
      % Assign to the subscripted beam

      if numel(subs) > 1
        if subs{1} == 1
          subs = subs(2:end);
        end
        assert(numel(subs) == 1, 'Only 1-D indexing supported for now');
      end

      assert(isempty(rem), 'Assignment to parts of beams not supported');
      if isempty(other)
        % Delete data
        beam.a(:, subs{:}) = [];
        beam.b(:, subs{:}) = [];

      else
        % Ensure we have a plane wave
        if ~isa(other, 'ott.beam.vswf.Bsc')
          other = ott.beam.vswf.Bsc(other);
        end

        % Ensure array sizes match
        beam.Nmax = max(beam.Nmax, other.Nmax);
        other.Nmax = beam.Nmax;

        beam.a(:, subs{:}) = other.a;
        beam.b(:, subs{:}) = other.b;
      end
    end

    function varargout = forcetorqueInternal(ibeam, other, varargin)
      % Calculate change in momentum between beams
      %
      % Usage
      %   [f, t, s] = ibeam.forcetorque(sbeam, ...) calculates the force,
      %   torque and spin between the incident beam ``ibeam`` and
      %   scattered beam ``sbeam``.
      %   Outputs 3xN matrix depending on the number of beams and
      %   other optional arguments, see bellow for more details.
      %
      %   [f, t, s] = beam.forcetorque(Tmatrix, ...) as above but first
      %   calculates the scattered beam using the given T-matrix.
      %
      %   [fx, fy, fz, tx, ty, tz, sx, sy, sz] = beam.forcetorque(...) as
      %   above but returns the x, y and z components as separate vectors.
      %
      % Optional named arguments
      %   - position (3xN numeric) -- Distance to translate beam before
      %     calculating the scattered beam using the T-matrix.
      %     Default: ``[]``.
      %   - rotation (3x3N numeric) -- Angle to rotate beam before
      %     calculating the scattered beam using the T-matrix.
      %     Inverse rotation is applied to scattered beam, effectively
      %     rotating the particle.
      %     Default: ``[]``.
      %
      % For details about position and rotation, see :meth:`scatter`.
      %
      % This uses mathematical result of Farsund et al., 1996, in the form of
      % Chricton and Marsden, 2000, and our standard T-matrix notation S.T.
      % E_{inc}=sum_{nm}(aM+bN);

      % Parse inputs
      [ibeam, sbeam, incN] = ibeam.forcetorqueParser(other, varargin{:});

      % Dispatch to other methods to calculate quantities
      force = ibeam.force(sbeam);
      torque = [];
      spin = [];
      if nargout > 1
        torque = ibeam.torque(sbeam);
        if nargout > 2
          spin = ibeam.spin(sbeam);
        end
      end

      % Combine incoherent beams
      if incN > 1
        force = squeeze(sum(reshape(force, 3, incN, []), 2));
        torque = squeeze(sum(reshape(torque, 3, incN, []), 2));
        spin = squeeze(sum(reshape(spin, 3, incN, []), 2));
      end

      % Package outputs
      if nargout <= 3
        varargout{1} = force;
        varargout{2} = torque;
        varargout{3} = spin;
      else
        varargout{1:3} = {force(1, :), force(2, :), force(3, :)};
        varargout{4:6} = {torque(1, :), torque(2, :), torque(3, :)};
        varargout{7:9} = {spin(1, :), spin(2, :), spin(3, :)};
      end
    end

    function varargout = forceInternal(ibeam, other, varargin)
      % Calculate change in linear momentum between beams.
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   force = ibeam.force(...)
      %   [fx, fy, fz] = ibeam.force(...)

      % Parse inputs
      [ibeam, sbeam, incN] = ibeam.forcetorqueParser(other, varargin{:});

      % Get the abpq terms for the calculation
      [a, b, p, q, n, m, ...
        anp1, bnp1, pnp1, qnp1, ...
        amp1, bmp1, pmp1, qmp1, ...
        anp1mp1, bnp1mp1, pnp1mp1, qnp1mp1, ...
        anp1mm1, bnp1mm1, pnp1mm1, qnp1mm1] = ...
      ibeam.find_abpq_force_terms(sbeam);

      % Calculate the Z force
      Az=m./n./(n+1).*imag(-(a).*conj(b)+conj(q).*(p));
      Bz=1./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ... %.*n
          .*imag(anp1.*conj(a)+bnp1.*conj(b)-(pnp1).*conj(p) ...
          -(qnp1).*conj(q));
      fz=2*sum(Az+Bz);

      % Calculate the XY force
      Axy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)) ...
          .*(conj(pmp1).*q - conj(amp1).*b - conj(qmp1).*p + conj(bmp1).*a);
      Bxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ... %sqrt(n.*)
          ( sqrt((n+m+1).*(n+m+2)) .* ( p.*conj(pnp1mp1) + q.* ...
          conj(qnp1mp1) -a.*conj(anp1mp1) -b.*conj(bnp1mp1)) + ...
          sqrt((n-m+1).*(n-m+2)) .* (pnp1mm1.*conj(p) + qnp1mm1.* ...
          conj(q) - anp1mm1.*conj(a) - bnp1mm1.*conj(b)) );

      fxy=sum(Axy+Bxy);
      fx=real(fxy);
      fy=imag(fxy);

      % Ensure things are full
      fx = full(fx);
      fy = full(fy);
      fz = full(fz);

      % Package output
      if nargout == 3
        varargout{1:3} = {fx, fy, fz};
      else
        varargout{1} = [fx(:) fy(:) fz(:)].';
      end
    end

    function varargout = torqueInternal(ibeam, other, varargin)
      % Calculate change in angular momentum between beams
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   torque = ibeam.torque(...)
      %   [tx, ty, tz] = ibeam.torque(...)

      % Parse inputs
      [ibeam, sbeam, incN] = ibeam.forcetorqueParser(other, varargin{:});

      % Get the abpq terms for the calculation
      [a, b, p, q, n, m, ~, ~, ~, ~, ...
        amp1, bmp1, pmp1, qmp1] = ...
      ibeam.find_abpq_force_terms(sbeam);

      tz=sum(m.*(a.*conj(a)+b.*conj(b)-p.*conj(p)-q.*conj(q)));

      txy=sum(sqrt((n-m).*(n+m+1)).*(a.*conj(amp1)+...
        b.*conj(bmp1)-p.*conj(pmp1)-q.*conj(qmp1)));
      tx=real(txy);
      ty=imag(txy);

      % Combine incoherent beams
      if incN > 1
        tx = sum(reshape(tx, incN, []), 1);
        ty = sum(reshape(ty, incN, []), 1);
        tz = sum(reshape(tz, incN, []), 1);
      end

      % Ensure things are full
      tx = full(tx);
      ty = full(ty);
      tz = full(tz);

      % Package output
      if nargout == 3
        varargout{1:3} = {tx, ty, tz};
      else
        varargout{1} = [tx(:) ty(:) tz(:)].';
      end
    end

    function varargout = spinInternal(ibeam, other, varargin)
      % Calculate change in spin between beams
      % For details on usage/arguments see :meth:`forcetorque`.
      %
      % Usage
      %   torque = ibeam.torque(...)
      %   [tx, ty, tz] = ibeam.torque(...)

      % Parse inputs
      [ibeam, sbeam, incN] = ibeam.forcetorqueParser(other, varargin{:});

      % Get the abpq terms for the calculation
      [a, b, p, q, n, m, ...
        anp1, bnp1, pnp1, qnp1, ...
        amp1, bmp1, pmp1, qmp1, ...
        anp1mp1, bnp1mp1, pnp1mp1, qnp1mp1, ...
        anp1mm1, bnp1mm1, pnp1mm1, qnp1mm1] = ...
      ibeam.find_abpq_force_terms(sbeam);

      Cz=m./n./(n+1).*(-(a).*conj(a)+conj(q).*(q)-(b).*conj(b)+conj(p).*(p));
      Dz=-2./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ...
            .*real(anp1.*conj(b)-bnp1.*conj(a)-(pnp1).*conj(q) ...
            +(qnp1).*conj(p));

      sz = sum(Cz+Dz);

      Cxy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)).* ...
          (conj(pmp1).*p - conj(amp1).*a + conj(qmp1).*q - conj(bmp1).*b);
      Dxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ...
            ( (sqrt((n+m+1).*(n+m+2)) .* ...
            ( p.*conj(qnp1mp1) - q.* conj(pnp1mp1) - ...
            a.*conj(bnp1mp1) +b.*conj(anp1mp1))) + ...
            (sqrt((n-m+1).*(n-m+2)) .* ...
            (pnp1mm1.*conj(q) - qnp1mm1.*conj(p) ...
            - anp1mm1.*conj(b) + bnp1mm1.*conj(a))) );

      sxy=sum(Cxy+Dxy);
      sy=real(sxy);
      sx=imag(sxy);

      % Ensure things are full
      sx = full(sx);
      sy = full(sy);
      sz = full(sz);

      % Combine incoherent beams
      if incN > 1
        sx = sum(reshape(sx, incN, []), 1);
        sy = sum(reshape(sy, incN, []), 1);
        sz = sum(reshape(sz, incN, []), 1);
      end

      % Package output
      if nargout == 3
        varargout{1:3} = {sx, sy, sz};
      else
        varargout{1} = [sx(:) sy(:) sz(:)].';
      end
    end

    function [E, data] = efieldInternal(beam, xyz, varargin)
      % Method used by efield(xyz)
      [E, ~, data] = beam.ehfield(xyz, varargin{:});
    end
    function [H, data] = hfieldInternal(beam, xyz, varargin)
      % Method used by hfield(xyz)
      [~, H, data] = beam.ehfield(xyz, varargin{:});
    end

    function [E, data] = efarfieldInternal(beam, rtp, varargin)
      % Method used by efield(xyz)
      %
      % Optional parameters
      %   - basis (enum) -- Basis to use for visualisation.  Can either
      %     be 'incoming' or 'outgoing'.  Defaults to the beam basis
      %     unless 'regular', in which case uses 'incoming' basis.
      %
      % See :class:`ott.beam.Beam` for usage/parameters.

      [beam, unmatched] = beam.farfieldVisHelper(varargin{:});
      [E, ~, data] = beam.ehfarfield(rtp, unmatched{:});
    end
    function [H, data] = hfarfieldInternal(beam, rtp, varargin)
      % Method used by hfield(xyz)
      %
      % Optional parameters
      %   - basis (enum) -- Basis to use for visualisation.  Can either
      %     be 'incoming' or 'outgoing'.  Defaults to the beam basis
      %     unless 'regular', in which case uses 'incoming' basis.
      %
      % See :class:`ott.beam.Beam` for usage/parameters.

      [beam, unmatched] = beam.farfieldVisHelper(varargin{:});
      [~, H, data] = beam.ehfarfield(rtp, unmatched{:});
    end

    function [a, b, p, q, n, m, ...
      anp1, bnp1, pnp1, qnp1, ...
      amp1, bmp1, pmp1, qmp1, ...
      anp1mp1, bnp1mp1, pnp1mp1, qnp1mp1, ...
      anp1mm1, bnp1mm1, pnp1mm1, qnp1mm1] = ...
    find_abpq_force_terms(ibeam, sbeam)
      % Find the terms required for force/torque calculations
      %
      % From forcetorque function in OTTv1

      % Ensure beams are the same size
      if ibeam.Nmax > sbeam.Nmax
        sbeam.Nmax = ibeam.Nmax;
      elseif ibeam.Nmax < sbeam.Nmax
        ibeam.Nmax = sbeam.Nmax;
      end

      % Ensure the beam is incoming-outgoing
      if isa(sbeam, 'ott.beam.abstract.Scattered')
        if ~strcmpi(sbeam.type, 'total')
          warning('Scattered beam type should be set to ''total''');
        end
      end

      % Get the relevent beam coefficients
      [a, b] = ibeam.getCoefficients();
      [p, q] = sbeam.getCoefficients();
      [n, m] = ibeam.getModeIndices();

      nmax=ibeam.Nmax;

      b=1i*b;
      q=1i*q;

      addv=zeros(2*nmax+3,1);

      at=[a;repmat(addv, 1, size(a, 2))];
      bt=[b;repmat(addv, 1, size(b, 2))];
      pt=[p;repmat(addv, 1, size(p, 2))];
      qt=[q;repmat(addv, 1, size(q, 2))];

      ci=ott.utils.combined_index(n,m);

      %these preserve order and number of entries!
      np1=2*n+2;
      cinp1=ci+np1;
      cinp1mp1=ci+np1+1;
      cinp1mm1=ci+np1-1;
      cimp1=ci+1;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %this is for m+1... if m+1>n then we'll ignore!
      kimp=(m>n-1);

      anp1=at(cinp1, :);
      bnp1=bt(cinp1, :);
      pnp1=pt(cinp1, :);
      qnp1=qt(cinp1, :);

      anp1mp1=at(cinp1mp1, :);
      bnp1mp1=bt(cinp1mp1, :);
      pnp1mp1=pt(cinp1mp1, :);
      qnp1mp1=qt(cinp1mp1, :);

      anp1mm1=at(cinp1mm1, :);
      bnp1mm1=bt(cinp1mm1, :);
      pnp1mm1=pt(cinp1mm1, :);
      qnp1mm1=qt(cinp1mm1, :);

      amp1=at(cimp1, :);
      bmp1=bt(cimp1, :);
      pmp1=pt(cimp1, :);
      qmp1=qt(cimp1, :);

      amp1(kimp, :)=0;
      bmp1(kimp, :)=0;
      pmp1(kimp, :)=0;
      qmp1(kimp, :)=0;

      a=a(ci, :);
      b=b(ci, :);
      p=p(ci, :);
      q=q(ci, :);

    end

    function [ibeam, sbeam, incN] = forcetorqueParser(ibeam, sbeam, varargin)
      % Input parser for forcetorque and related methods

      assert(isa(ibeam, 'ott.beam.vswf.Bsc'), 'ibeam must be Bsc');
      assert(isa(sbeam, 'ott.beam.vswf.Bsc'), 'sbeam must be Bsc');

      % Combine beams if coherent
      if strcmpi(ibeam.array_type, 'coherent')
        ibeam = sum(ibeam);
      end
      if strcmpi(sbeam.array_type, 'coherent')
        sbeam = sum(sbeam);
      end

      % Ensure beams are total beams
      if isa(ibeam, 'ott.beam.abstract.Scattered')
        ibeam = ibeam.total_beam;
      end
      if isa(sbeam, 'ott.beam.abstract.Scattered')
        sbeam = sbeam.total_beam;
      end

      % Get the number of beams (for combination)
      assert(numel(ibeam) == 1 || numel(sbeam) == 1 ...
          || numel(ibeam) == numel(sbeam), ...
          'Number of incident and scattered beams must be 1 or matching');
      Nbeams = max(numel(ibeam), numel(sbeam));

      % Handle the combine argument
      % TODO: This should use the array functions
      incN = 1;
      if strcmpi(ibeam.array_type, 'incoherent')
        % Return how many terms need to be combined later
        incN = Nbeams;
      end
    end

    function [beam, unmatched] = farfieldVisHelper(beam, varargin)
      % Helper for far-field visualisation functions

      p = inputParser;
      p.addParameter('basis', []);
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      unmatched = ott.utils.unmatchedArgs(p);

      % Change basis
      if isempty(p.Results.basis) && strcmpi(beam.basis, 'regular')
        beam.basis = 'incoming';
      elseif ~isempty(p.Results.basis)
        beam.basis = p.Results.basis;
      end
    end
  end

  methods % Getters/setters
    % power
    % Nmax
    % basis
    % absdz

    function nmax = get.Nmax(beam)
      % Calculates Nmax from the current size of the beam coefficients
      nmax = ott.utils.combined_index(size(beam.a, 1));
    end
    function beam = set.Nmax(beam, nmax)
      % Resizes the beam vectors (a,b)
      beam = beam.setNmax(nmax);
    end

    function p = get.power(beam)
      % get.power calculate the power of the beam
      % TODO: Shouldn't this change for different array types
      p = full(sum(abs(beam.a).^2 + abs(beam.b).^2));
    end
    function beam = set.power(beam, val)
      % set.power set the beam power
      beam = sqrt(val / beam.power) * beam;
    end

    function beam = set.basis(beam, val)
      assert(any(strcmpi(val, {'incoming', 'outgoing', 'regular'})), ...
        'ott:beam:vswf:Bsc:set_basis:invalid_value', ...
        'basis must be one of ''incoming'' ''outgoing'' or ''regular''');

      beam.basis = val;
    end

    function beam = set.absdz(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'absdz must be numeric scalar');
      beam.absdz = val;
    end
  end
end
