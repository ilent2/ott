classdef Bsc < ott.optics.beam.Beam
% Class representing beam shape coefficients.
% Inherits from :class:`ott.optics.beam.Beam`.
%
% Any units can be used for the properties as long as they are
% consistent in all specified properties.  Calculated quantities
% will have these units.
%
% Properties
%   - a           --  Beam shape coefficients a vector
%   - b           --  Beam shape coefficients b vector
%   - basis       --  VSWF beam basis (incoming, outgoing or regular)
%   - omega       --  Angular frequency of beam [2*pi/T]
%   - wavelength  --  Wavelength of beam [L]
%   - dz          --  Absolute cumulative distance the beam has moved
%
% Dependent properties
%   - Nmax        --  Truncation number for VSWF coefficients
%   - power       --  Power of the beam [M*L^2/S^3]
%   - Nbeams      --  Number of beams in this Bsc object
%   - wavenumber  --  Wavenumber in the medium [2*pi/L]
%   - speed       --  Speed of beam in medium [L/T]
%
% Methods
%   - append      --  Joins two beam objects together
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
%
% Static methods:
%   - make_beam_vector -- Convert output of bsc_* functions to beam coefficients
%
% See also Bsc, ott.BscPmGauss, ott.BscPlane.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    a           % Beam shape coefficients a vector
    b           % Beam shape coefficients b vector

    dz          % Absolute cumulative distance the beam has moved

    % These can't be tracked using Matrix translation/rotations
    %offset      % Offset applied to beam using translate functions
    %direction   % Direction of beam applied using rotation functions
  end

  properties
    basis       % VSWF beam basis (incoming, outgoing or regular)
    omega       % Angular frequency of beam
  end

  properties (Dependent)
    Nmax        % Truncation number for VSWF coefficients
    Nbeams      % Number of beams in this Bsc object

    speed       % Speed of beam in medium
  end

  methods (Static)
    function [a, b, n, m] = make_beam_vector(a, b, n, m, Nmax)
      %MAKE_BEAM_VECTOR converts output of bsc_* functions to sparse vectors

      % Check we have some modes
      if isempty(n)
        warning('OTT:BSC:make_beam_vector:no_modes', ...
            'No modes in beam or all zero modes.');
      end

      % Handle default value for Nmax
      if nargin < 5
        if isempty(n)
          Nmax = 0;
        else
          Nmax = max(n);
        end
      end

      total_orders = ott.utils.combined_index(Nmax, Nmax);
      ci = ott.utils.combined_index(n, m);
      nbeams = size(a, 2);

      [ci, cinbeams] = ndgrid(ci, 1:nbeams);

      a = sparse(ci, cinbeams, a, total_orders, nbeams);
      b = sparse(ci, cinbeams, b, total_orders, nbeams);

      [n, m] = ott.utils.combined_index(1:Nmax^2+2*Nmax);
      n = n.';
      m = m.';
    end

    function k_medium = parser_k_medium(p, default)
      %PARSER_K_MEDIUM helper to get k_medium from a parser object

      if ~isempty(p.Results.k_medium)
        k_medium = p.Results.k_medium;
      elseif ~isempty(p.Results.wavelength_medium)
        k_medium = 2.0*pi/p.Results.wavelength_medium;
      elseif ~isempty(p.Results.index_medium)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_medium');
        end
        k_medium = p.Results.index_medium*2.0*pi/p.Results.wavelength0;
      elseif nargin == 2
        k_medium = default;
      else
        error('Unable to determine k_medium from inputs');
      end
    end
  end

  methods (Access=protected)
    
    function [A, B] = translateZ_type_helper(beam, z, Nmax)
      % Determine the correct type to use in ott.utils.translate_z
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
    
  end

  methods
    function beam = Bsc(varargin)
      % Construct a new beam object
      %
      % Usage
      %   beam = Bsc(a, b, basis, ...) constructs a new beam vector.
      %   Useful if you have a specific set of a/b coefficients that you
      %   want to wrap in a beam object.
      %
      %   beam = Bsc(a, b)  As above, but uses 'regular' as the basis.
      %
      %   beam = Bsc()  Constructs an empty beam object.
      %   beam = Bsc(bsc)  Construct a copy of the existing bsc object.
      %
      % Parameters
      %   - a,b (numeric) -- Vectors of VSWF coefficients
      %   - bsc (vswf.bsc.Bsc) -- An existing BSC object
      %   - basis (enum) -- VSWF basis: incoming, outgoing or regular
      %
      % Optional named arguments
      %   - k_medium  n -- Wavenumber in medium (default: 2*pi)
      %   - omega     n -- Angular frequency (default: 2*pi)
      %   - dz        n -- Initial displacement of the beam (default: 0)
      %   - like   beam -- Construct this beam to be like another beam

      p = inputParser;
      p.addOptional('a', [], ...
        @(x) isnumeric(x) || isa(x, 'ott.optics.vswf.bsc.Bsc'));
      p.addOptional('b', [], @isnumeric);
      p.addOptional('basis', 'regular', ...
        @(x) any(strcmpi(x, {'incoming', 'outgoing', 'regular'})));
      p.addParameter('like', []);
      p.addParameter('k_medium', 2.0*pi);
      p.addParameter('omega', 2*pi);
      p.addParameter('dz', 0.0);
      p.parse(varargin{:});

      a = p.Results.a;
      b = p.Results.b;
      like = p.Results.like;

      % Copy all data from an existing BSC
      if isa(a, 'ott.optics.vswf.bsc.Bsc')
        assert(isempty(b), 'Too many arguments supplied');

        like = a;
        a = like.a;
        b = like.b;
      end

      % Copy all beam data (except coefficients) from like
      if ~isempty(like)
        beam.dz = like.dz;
        beam.wavenumber = like.wavenumber;
        beam.omega = like.omega;
        beam.basis = like.basis;
      else
        beam.dz = p.Results.dz;
        beam.wavenumber = p.Results.k_medium;
        beam.omega = p.Results.omega;
        beam.basis = p.Results.basis;
      end

      % Check size of a and b and assign
      assert(all(size(a) == size(b)), 'size of a and b must match');
      assert(isempty(a) || ...
          (size(a, 1) >= 3 && ...
          sqrt(size(a, 1)+1) == floor(sqrt(size(a, 1)+1))), ...
        'number of multipole terms must be empty, 3, 8, 15, 24, ...');
      beam.a = a;
      beam.b = b;
    end

    function beam = append(beam, other)
      % APPEND joins two beam objects together
      %
      % Usaage
      %   beam = beam.append(other)
      %   beam = beam.append([beam1, beam2, ...])
      
      % Append array of beams
      if numel(other) > 1
        for ii = 1:numel(other)
          beam = beam.append(other(ii));
        end
        return;
      end

      % Append single beam
      if beam.Nbeams == 0
        % Copy the other beam, preserves properties
        beam = other;
      else
        beam.Nmax = max(beam.Nmax, other.Nmax);
        other.Nmax = beam.Nmax;
        beam.a = [beam.a, other.a];
        beam.b = [beam.b, other.b];
      end
    end
    
    function varargout = paraxialFarfield(beam, varargin)
      % Calcualtes fields in the paraxial far-field
      %
      % Usage
      %   [E, H] = beam.paraxialFarfield(...)
      %
      % Optional named arguments
      %   - 'mapping' (enum) -- mapping to paraxial far-field
      %   - 'calcE'   bool   calculate E field (default: true)
      %   - 'calcH'   bool   calculate H field (default: nargout == 2)
      %   - 'saveData' bool  save data for repeated calculation (default: false)
      %   - 'data'    data   data saved for repeated calculation.
      
      p = inputParser;
      p.addParameter('calcE', true);
      p.addParameter('calcH', nargout >= 2);
      p.addParameter('saveData', false);
      p.addParameter('data', []);
      p.addParameter('mapping', 'sin');
      p.addParameter('size', [50, 50]);
      p.addParameter('thetaMax', pi/2);
      p.addParameter('direction', 'pos');
      p.parse(varargin{:});
      
      % Calculate image locations
      xrange = linspace(-1, 1, p.Results.size(1));
      yrange = linspace(-1, 1, p.Results.size(2));
      [xx, yy] = meshgrid(xrange, yrange);

      % Calculate spherical coordinates for pixels
      phi = atan2(yy, xx);
      rr = sqrt(xx.^2 + yy.^2);
      switch p.Results.mapping
        case 'sin'
          theta = asin(rr);
        case 'tan'
          theta = atan(rr);
        case 'theta'
          theta = rr;
        otherwise
          error('Unknown mapping argument value, must be sin, tan or theta');
      end
      
      % Only include points within NA range
      thetaMax = Inf;
      if ~isempty(p.Results.thetaMax)
        thetaMax = p.Results.thetaMax;
      end

      % Determine if the points need calculating
      pinside = imag(theta) == 0 & theta < thetaMax;
      iphi = phi(pinside);
      itheta = theta(pinside);
      
      if strcmpi(p.Results.direction, 'neg')
        itheta = pi - itheta;
      elseif ~strcmpi(p.Results.direction, 'pos')
        error('Direction must be ''pos'' or ''neg''');
      end

      % Calculate the electric field in the farfield
      [E, H, data] = beam.farfield(itheta(:), iphi(:), ...
        'saveData', p.Results.saveData, 'data', p.Results.data, ...
        'calcE', p.Results.calcE, 'calcH', p.Results.calcH);
      
      if nargout >= 1
        
        if p.Results.calcE
          % Generate the requested field
          dEt = beam.GetVisualisationData('Et', [], ...
            [itheta, iphi, ones(size(iphi))], [], E.');
          dEp = beam.GetVisualisationData('Ep', [], ...
            [itheta, iphi, ones(size(iphi))], [], E.');

          Et = zeros(p.Results.size);
          Et(pinside) = dEt;
          Ep = zeros(p.Results.size);
          Ep(pinside) = dEp;

          varargout{1} = Et;
          varargout{1}(:, :, 2) = Ep;
        end
        
        if nargout >= 2 && p.Results.calcH
          % Generate the requested field
          dHt = beam.GetVisualisationData('Et', [], ...
            [itheta, iphi, ones(size(iphi))], [], H.');
          dHp = beam.GetVisualisationData('Ep', [], ...
            [itheta, iphi, ones(size(iphi))], [], H.');

          Ht = zeros(p.Results.size);
          Ht(pinside) = dHt;
          Hp = zeros(p.Results.size);
          Hp(pinside) = dHp;

          varargout{2} = Ht;
          varargout{2}(:, :, 2) = Hp;
        end
      end
      
      if nargout >= 3
        varargout{3} = data;
      end
    end

    function [E, H, data] = farfield(beam, theta, phi, varargin)
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

    function [E, H, data] = emFieldRtp(beam, rtp, varargin)
      % Calculates the E and H field at specified locations
      %
      % [E, H] = beam.emFieldRtp(rtp, ...) calculates the complex field
      % at locations xyz (3xN matrix of spherical coordinates).
      % Returns 3xN matrices for the E and H field at these locations.
      %
      % Optional named arguments:
      %    'calcE'   bool   calculate E field (default: true)
      %    'calcH'   bool   calculate H field (default: nargout == 2)
      %    'saveData' bool  save data for repeated calculation (default: false)
      %    'data'    data   data saved for repeated calculation.
      %    'coord'   str    coordinates to use for calculated field
      %       'cartesian' (default) or 'spherical'
      %
      % If either calcH or calcE is false, the function still returns
      % E and H as matrices of all zeros for the corresponding field.

      p = inputParser;
      p.addParameter('calcE', true);
      p.addParameter('calcH', nargout >= 2);
      p.addParameter('saveData', false);
      p.addParameter('data', []);
      p.parse(varargin{:});

      % Scale the locations by the wave number (unitless coordinates)
      rtp(1, :) = rtp(1, :) * abs(beam.wavenumber);

      % Get the indices required for the calculation
      [n,m]=ott.utils.combined_index(find(abs(beam.a)|abs(beam.b)));
      nm = [ n; m ];

      ci = ott.utils.combined_index(n, m);
      [a, b] = beam.getCoefficients(ci);

      % Calculate the fields
      % TODO: Should this function be moved to bsc.Static or bsc.private?
      % TODO: Should this function take 3xN instead of Nx3
      [E, H, data] = ott.utils.emField(rtp.', beam.basis, nm, [a; b], ...
          'saveData', p.Results.saveData, ...
          'data', p.Results.data, ...
          'calcE', p.Results.calcE, 'calcH', p.Results.calcH);
      E = E.';
      H = H.';

      % Package output
      E = ott.utils.FieldVector(rtp, E, 'spherical');
      H = ott.utils.FieldVector(rtp, H, 'spherical');
    end

    function [E, H] = ehfield(beam, xyz)
      % Calculate E and H field
      %
      % Usage
      %   [E, H] = beam.ehfield(xyz)
      %   Calculates the fields at the specified locations (3xN matrix).
      [E, H, data] = beam.emFieldXyz(xyz);
    end

    function [E, H, data] = emFieldXyz(beam, xyz, varargin)
      %EMFIELDXYZ calculates the E and H field at specified locations
      %
      % [E, H] = beam.emFieldXyz(xyz, ...) calculates the complex field
      % at locations xyz (3xN matrix of Cartesian coordinates).
      % Returns 3xN matrices for the E and H field at these locations.
      %
      % Optional named arguments:
      %    'calcE'   bool   calculate E field (default: true)
      %    'calcH'   bool   calculate H field (default: nargout == 2)
      %    'saveData' bool  save data for repeated calculation (default: false)
      %    'data'    data   data saved for repeated calculation.
      %    'coord'   str    coordinates to use for calculated field
      %       'cartesian' (default) or 'spherical'
      %
      % If either calcH or calcE is false, the function still returns
      % E and H as matrices of all zeros for the corresponding field.

      rtp = ott.utils.xyz2rtp(xyz);
      [E, H, data] = beam.emFieldRtp(rtp, varargin{:});
    end
    
    function varargout = visualiseFarfieldSlice(beam, phi, varargin)
      % Generate a 2-D scattering plot of the far-field
      %
      % beam.visualiseFarfieldSlice(phi) display a visualisation
      % of the farfield.
      %
      % [theta, I] = beam.visualiseFarfieldSlice(phi) calculate data
      % for the visualisation.  Pass showVisualisation, true to show.
      
      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('normalise', false);
      p.addParameter('ntheta', 100);
      p.addParameter('showVisualisation', nargout == 0);
      p.parse(varargin{:});
      
      ptheta = linspace(0, 2*pi, p.Results.ntheta);
      
      % TODO: Other field types

      % Calculate electric field
      [E, ~] = beam.farfield(ptheta, phi);
      
      % Calculate desired field
      [rtp{1:3}] = ott.utils.matchsize(0, ptheta(:), phi);
      I = beam.GetVisualisationData(p.Results.field, [], [rtp{1}, rtp{2}, rtp{3}], [], E.');
%       I = sum(abs(E).^2, 1);
      
      if p.Results.normalise
        I = I ./ max(abs(I(:)));
      end

      % Setup outputs
      if nargout == 2
        varargout{1} = theta;
        varargout{2} = I;
      end
      
      % Display visualisation
      if p.Results.showVisualisation
        polarplot(ptheta, I);
      end
      
    end
    
    function visualiseFarfieldSphere(beam, varargin)
      % Generate a spherical surface visualisation of the far-field
      %
      % beam.visualiseFarfieldSphere(phi)
      %
      % Optional named arguments:
      %   npts      num   Number of points to use for sphere surface
      %   normalise bool  If intensity values should be normalised to 1
      %   type      str   Type of visualisation to produce.
      %       sphere    (default) draw a sphere with intensity as color
      %       3dpolar   scale the radius by the intensity
      
      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('npts', 100);
      p.addParameter('normalise', false);
      p.addParameter('type', 'sphere');
      p.parse(varargin{:});
      
      % build grid:
      [x,y,z]=sphere(p.Results.npts);

      % generate angular points for farfield:
      [~,theta,phi]=ott.utils.xyz2rtp(x,y,z);

      % find far-field in theta, phi:
      [E,~]=beam.farfield(theta(:),phi(:));
      
      % Calculate the requested field
      dataout = beam.GetVisualisationData(p.Results.field, [], ...
        [theta(:), phi(:), ones(size(phi(:)))], [], E.');

      % Reshape to match the input
      I=reshape(dataout,size(x));
      
      if p.Results.normalise
        I = I ./ max(abs(I(:)));
      end
      
      switch p.Results.type
        case 'sphere'
          surf(x,y,z,I,'facecolor','interp','edgecolor','none');
        case '3dpolar'
          surf(abs(I).*x,abs(I).*y,abs(I).*z,I,...
              'facecolor','interp','edgecolor','none');
        otherwise
          error('Unknown visualisation type');
      end

      zlabel('Z');
      xlabel('X');
      ylabel('Y');
      view(50, 20);
      axis equal;
    end

    function varargout = visualiseFarfield(beam, varargin)
      % Create a 2-D visualisation of the farfield of the beam
      %
      % visualiseFarfield(...) displays an image of the farfield in
      % the current figure window.
      %
      % im = visualiseFarfield(...) returns a 2-D image of the farfield.
      %
      % [im, data] = visualiseFarfield(..., 'saveData', true) returns the
      % saved data that can be used for repeated calculation.
      %
      %    TODO: Should the data object instead be a callable object?
      %     This would make the interface simpler.
      %
      % Optional named arguments:
      %     'size'    [ x, y ]    Size of the image
      %     'direction'  dir      Hemisphere string ('pos' or 'neg'),
      %        2-vector (theta, phi) or 3x3 rotation matrix.
      %     'field'   type        Type of field to calculate
      %     'mapping' map         Mapping from sphere to plane ('sin', 'tan')
      %     'range'   [ x, y ]    Range of points to visualise
      %    'saveData' bool  save data for repeated calculation (default: false)
      %    'data'    data   data saved for repeated calculation.
      %    'thetaMax' num   maximum theta angle to include in image
      %    'showVisualisation'  bool   show the visualisation in the
      %       current figure (default: nargout == 0).

      p = inputParser;
      p.addParameter('size', [80, 80]);
      p.addParameter('direction', 'pos');
      p.addParameter('field', 'irradiance');
      p.addParameter('mapping', 'sin');
      p.addParameter('range', [1, 1]);
      p.addParameter('saveData', nargout == 2);
      p.addParameter('data', []);
      p.addParameter('thetaMax', []);
      p.addParameter('showVisualisation', nargout == 0);
      p.parse(varargin{:});
      
      % If direction is a vector, rotate to that direction
      if ~ischar(p.Results.direction)
        dir = p.Results.direction;
        if numel(dir) == 2
          rbeam = beam.rotateYz(dir(1), dir(2));
        elseif all(size(dir) == [3, 3])
          rbeam = beam.rotate(dir);
        else
          error('OTT:BSC:visualiseFarfield:bad_direction', ...
            'Direction must be char array or 2 element vector or 3x3 matrix');
        end
        
        [varargout{1:nargout}] = rbeam.visualiseFarfield(...
          'size', p.Results.size, 'direction', 'pos', ...
          'field', p.Results.field, 'mapping', p.Results.mapping, ...
          'range', p.Results.range, 'saveData', p.Results.saveData, ...
          'data', p.Results.data, 'thetaMax', p.Results.thetaMax, ...
          'showVisualisation', p.Results.showVisualisation);
        return;  % All done
      end

      % Calculate image locations
      xrange = linspace(-1, 1, p.Results.size(1))*p.Results.range(1);
      yrange = linspace(-1, 1, p.Results.size(2))*p.Results.range(2);
      [xx, yy] = meshgrid(xrange, yrange);

      % Calculate spherical coordinates for pixels
      phi = atan2(yy, xx);
      rr = sqrt(xx.^2 + yy.^2);
      switch p.Results.mapping
        case 'sin'
          theta = asin(rr);
        case 'tan'
          theta = atan(rr);
        case 'theta'
          theta = rr;
        otherwise
          error('Unknown mapping argument value, must be sin, tan or theta');
      end
      
      % Only include points within NA range
      thetaMax = Inf;
      if ~isempty(p.Results.thetaMax)
        thetaMax = p.Results.thetaMax;
      end

      % Determine if the points need calculating
      pinside = imag(theta) == 0 & theta < thetaMax;
      iphi = phi(pinside);
      itheta = theta(pinside);
      
      if strcmpi(p.Results.direction, 'neg')
        itheta = pi - itheta;
      elseif ~strcmpi(p.Results.direction, 'pos')
        error('Direction must be ''pos'' or ''neg''');
      end

      % Calculate the electric field in the farfield
      [ioutputE, ~, data] = beam.farfield(itheta(:), iphi(:), ...
        'saveData', p.Results.saveData, 'data', p.Results.data, ...
        'calcE', true, 'calcH', false);
      
      % Generate the requested field
      
      dataout = beam.VisualisationData(p.Results.field, ioutputE);

      % Pack the result into the images
      imout = zeros(p.Results.size);
      imout(pinside) = dataout;

      % Display the visualisation
      if p.Results.showVisualisation
        
        % Check the field is real
        if ~isreal(imout)
          error(['Unsupported field type for visualisation: ' p.Results.field]);
        end
        
        imagesc(xrange, yrange, imout);
        caxis([min(dataout), max(dataout)]);
        xlabel('X');
        ylabel('Y');
        axis image;
      end
      
      % Handle outputs
      if nargout == 1
        varargout{1} = imout;
      elseif nargout == 2
        varargout{1} = imout;
        varargout{2} = data;
      end
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
      
      assert(numel(beam) == 1, ...
        'Only single beam inputs supported for Bsc.visualise');

      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('axis', 'z');
      p.addParameter('combine', []);
      p.addParameter('range', ...
          [1,1].*ott.utils.nmax2ka(beam.Nmax)./abs(beam.wavenumber));
      p.addParameter('size', []);
      p.addParameter('plot_axes', []);
      p.parse(varargin{:});

      assert(isempty(p.Results.combine) || ...
          any(strcmpi(p.Results.combine, {'coherent', 'incoherent'})), ...
          'combine must be one of empty, ''coherent'' or ''incoherent''');

      % Get unmatched arguments
      unmatched = [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'];

      % Combine coherent beams
      if strcmpi(p.Results.combine, 'coherent') || beam.Nbeams == 1
        beam = sum(beam);
        [varargout{1:nargout}] = visualise@ott.optics.beam.Beam(...
            beam, unmatched{:}, 'range', p.Results.range, ...
            'size', p.Results.size, 'axis', p.Results.axis, ...
            'plot_axes', p.Results.plot_axes);
      else

        % Separate beams into beam array
        beam_array = beam.beam(1);
        for ii = 2:beam.Nbeams
          beam_array(ii) = beam.beam(ii);
        end

        % Generate the data
        imout = visualise@ott.optics.beam.Beam(...
            beam_array, unmatched{:}, 'range', p.Results.range, ...
            'size', p.Results.size, 'axis', p.Results.axis, ...
            'plot_axes', []);

        % Combine incoherently
        if strcmpi(p.Results.combine, 'incoherent')
          imout = sum(im, 3);
        end

        % Get data for plot
        default_sz = [80, 80];
        [xrange, yrange, ~] = beam.visualiseGetRange(p.Results.range, ...
            p.Results.size, default_sz);
        [~, labels] = beam.visualiseGetXyz([], [], [], p.Results.axis);

        % Call visualisation helper
        beam.visualiseShowPlot(nargout, p.Results.plot_axes, imout, ...
            {xrange, yrange}, labels);

        % Assign outputs if requested
        if nargout == 1
          varargout{1} = imout;
        end

      end
    end

    function speed = get.speed(beam)
      % Get the speed of the beam in medium
      speed = beam.omega / beam.wavenumber;
    end

    function beam = set.basis(beam, basis)
      % Set the beam type, checking it is a valid type first
      if ~any(strcmpi(basis, {'incoming', 'outgoing', 'regular'}))
        error('OTT:Bsc:set_basis:invalid_value', 'Invalid beam basis');
      end
      beam.basis = basis;
    end

    function nbeams = get.Nbeams(beam)
      % get.beams get the number of beams in this object
      nbeams = size(beam.a, 2);
    end

    function bsc = beam(bsc, idx)
      % BEAM get beams from a beam array object
      %
      % BEAM(idx) idx can be a linear index or a logical array.
      bsc.a = bsc.a(:, idx);
      bsc.b = bsc.b(:, idx);
    end

    function nmax = get.Nmax(beam)
      %get.Nmax calculates Nmax from the current size of the beam coefficients
      nmax = ott.utils.combined_index(size(beam.a, 1));
    end

    function beam = set.Nmax(beam, nmax)
      %set.Nmax resizes the beam vectors
      beam = beam.set_Nmax(nmax);
    end

    function nbeam = shrink_Nmax(beam, varargin)
      % SHRINK_NMAX reduces the size of the beam while preserving power

      p = inputParser;
      p.addParameter('tolerance', 1.0e-6);
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

    function beam = set_Nmax(beam, nmax, varargin)
      % SET_NMAX resize the beam, with additional options
      %
      % SET_NMAX(nmax) sets the beam nmax.
      %
      % SET_NMAX(..., 'tolerance', tol) use tol as the warning error
      % level tolerance for resizing the beam.
      %
      % SET_NMAX(..., 'powerloss', mode) action to take if a power
      % loss is detected.  Can be 'ignore', 'warn' or 'error'.

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
              warning('ott:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            elseif strcmpi(p.Results.powerloss, 'error')
              error('ott:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            else
              error('ott:Bsc:setNmax:truncation', ...
                'powerloss should be one of ignore, warn or error');
            end
          end
        end
      elseif size(beam.a, 1) < total_orders
        [arow_index,acol_index,aa] = find(beam.a);
        [brow_index,bcol_index,ba] = find(beam.b);
        beam.a = sparse(arow_index,acol_index,aa,total_orders,beam.Nbeams);
        beam.b = sparse(brow_index,bcol_index,ba,total_orders,beam.Nbeams);
      end
    end

    function beam = translate(beam, A, B)
      % TRANSLATE apply a translation using given translation matrices.
      %
      % TRANSLATE(A, B) applies the translation given by A, B.
      beam = [ A B ; B A ] * beam;
    end

    function [beam, A, B] = translateZ(beam, varargin)
      % Translate a beam along the z-axis.
      %
      % Units for the coordinates should be consistent with the
      % beam wave number (i.e., if the beam was created by specifying
      % wavelength in units of meters, distances here should also be
      % in units of meters).
      %
      % Usage
      %   tbeam = beam.translateZ(z) translates by a distance ``z``
      %   along the z axis.
      %
      %   [tbeam, A, B] = beam.translateZ(z) returns the translation matrices
      %   and the translated beam.  See also :meth:`+ott.Bsc.translate`.
      %
      %   [tbeam, AB] = beam.translateZ(z) returns the ``A, B`` matrices
      %   packed so they can be directly applied to a
      %   beam: ``tbeam = AB * beam``.
      %
      %   [...] = beam.translateZ(..., 'Nmax', Nmax) specifies the output
      %   beam ``Nmax``.  Takes advantage of not needing to calculate
      %   a full translation matrix.

      p = inputParser;
      p.addOptional('z', []);
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      if nargout ~= 1 && numel(p.Results.z) > 1
        error('Multiple output with multiple translations not supported');
      end

      if ~isempty(p.Results.z)
        z = p.Results.z;

        % Add a warning when the beam is translated outside nmax2ka(Nmax) 
        % The first time may be OK, the second time does not have enough
        % information.
        if any(beam.dz > ott.utils.nmax2ka(beam.Nmax)./beam.wavenumber)
          warning('ott:Bsc:translateZ:outside_nmax', ...
              'Repeated translation of beam outside Nmax region');
        end
        beam.dz = beam.dz + abs(z);

        % Convert to beam units
        z = z * beam.wavenumber / 2 / pi;

        ibeam = beam;
        beam = ott.optics.vswf.bsc.Bsc();

        for ii = 1:numel(z)
          [A, B] = ibeam.translateZ_type_helper(z(ii), [p.Results.Nmax, ibeam.Nmax]);
          beam = beam.append(ibeam.translate(A, B));
          beam.basis = 'regular';
        end
      else
        error('Wrong number of arguments');
      end

      % Pack the rotated matricies into a single ABBA object
      if nargout == 2
        A = [ A B ; B A ];
      end
    end

    function varargout = translateXyz(beam, varargin)
      % Translate the beam given Cartesian coordinates.
      %
      % Units for the coordinates should be consistent with the
      % beam wave number (i.e., if the beam was created by specifying
      % wavelength in units of meters, distances here should also be
      % in units of meters).
      %
      % Usage
      %   tbeam = beam.translateXyz(xyz) translate the beam to locations
      %   given by the ``xyz`` coordinates, where ``xyz`` is a 3xN matrix
      %   of coordinates.
      %
      %   tbeam = beam.translateXyz(Az, Bz, D)
      %   Translate the beam using z-translation and rotation matrices.
      %
      %   [tbeam, Az, Bz, D] = beam.translateXyz(...) returns the
      %   z-translation matrices ``Az, Bz``, the rotation matrix ``D``,
      %   and the translated beam ``tbeam``.
      %
      %   [tbeam, A, B] = beam.translateXyz(...) returns the translation
      %   matrices ``A, B`` and the translated beam.
      %
      %   [tbeam, AB] = beam.translateXyz(...) returns the translation
      %   matrices ``A, B`` packaged so they can be directly applied
      %   to a beam using ``tbeam = AB * beam``.
      %
      %   tbeam = beam.translateXyz(..., 'Nmax', Nmax) specifies the
      %   output beam ``Nmax``.  Takes advantage of not needing to
      %   calculate a full translation matrix.

      p = inputParser;
      p.addOptional('opt1', []);    % xyz or Az
      p.addOptional('opt2', []);    % [] or Bz
      p.addOptional('opt3', []);    % [] or D
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      if ~isempty(p.Results.opt1) && isempty(p.Results.opt2) ...
          && isempty(p.Results.opt3)
        xyz = p.Results.opt1;
        rtp = ott.utils.xyz2rtp(xyz);
        [varargout{1:nargout}] = beam.translateRtp(rtp, ...
            'Nmax', p.Results.Nmax);
      else
        [varargout{1:nargout}] = beam.translateRtp(varargin{:});
      end
    end

    function [beam, A, B, D] = translateRtp(beam, varargin)
      %TRANSLATERTP translate the beam given spherical coordinates
      %
      % TRANSLATERTP(rtp) translate the beam to locations given by
      % the xyz coordinates, where rtp is a 3xN matrix of coordinates.
      %
      % TRANSLATERTP(Az, Bz, D) translate the beam using
      % z-translation and rotation matricies.
      %
      % [beam, Az, Bz, D] = TRANSLATERTP(...) returns the z-translation
      % matrices, the rotation matrix D, and the translated beam.
      %
      % [beam, A, B] = TRANSLATERTP(...) returns the translation matrices
      % and the translated beam.
      %
      % [beam, AB] = TRANSLATERTP(...) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam.
      %
      % TRANSLATERTP(..., 'Nmax', Nmax) specifies the output beam Nmax.
      % Takes advantage of not needing to calculate a full translation matrix.

      p = inputParser;
      p.addOptional('opt1', []);    % rtp or Az
      p.addOptional('opt2', []);    % [] or Bz
      p.addOptional('opt3', []);    % [] or D
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      % Convert Nmax to a single number
      if numel(p.Results.Nmax) == 1
        oNmax = p.Results.Nmax;
      elseif numel(p.Results.Nmax) == 2
        oNmax = p.Results.Nmax(2);
      else
        error('Nmax must be 2 element vector or scalar');
      end

      % Handle input arguments
      if ~isempty(p.Results.opt1) && isempty(p.Results.opt2) ...
          && isempty(p.Results.opt3)

        % Assume first argument is rtp coordinates
        r = p.Results.opt1(1, :);
        theta = p.Results.opt1(2, :);
        phi = p.Results.opt1(3, :);

      elseif ~isempty(p.Results.opt1) && ~isempty(p.Results.opt2) ...
          && ~isempty(p.Results.opt3)

        % Rotation/translation is already computed, apply it
        A = p.Results.opt1;
        B = p.Results.opt2;
        D = p.Results.opt3;
        beam = beam.rotate('wigner', D);
        beam = beam.translate(A, B);
        beam = beam.rotate('wigner', D');
        return;
      else
        error('Not enough input arguments');
      end

      if numel(r) ~= 1 && nargout ~= 1
        error('Multiple output with multiple translations not supported');
      end

      % Only do the rotation if we need it
      if any((theta ~= 0 & abs(theta) ~= pi) | phi ~= 0)

        ibeam = beam;
        beam = ott.optics.vswf.bsc.Bsc();

        for ii = 1:numel(r)
          [tbeam, D] = ibeam.rotateYz(theta(ii), phi(ii), ...
              'Nmax', max(oNmax, ibeam.Nmax));
          [tbeam, A, B] = tbeam.translateZ(r(ii), 'Nmax', oNmax);
          beam = beam.append(tbeam.rotate('wigner', D'));
        end
      else
        dnmax = max(oNmax, beam.Nmax);
        D = speye(ott.utils.combined_index(dnmax, dnmax));

        % Replace rotations by 180 with negative translations
        idx = abs(theta) == pi;
        r(idx) = -r(idx);

        if numel(r) == 1
          [beam, A, B] = beam.translateZ(r, 'Nmax', oNmax);
        else
          beam = beam.translateZ(r, 'Nmax', oNmax);
        end
      end

      % Rotate the translation matricies
      if nargout == 3 || nargout == 2

        % The beam might change size, so readjust D to match
        sz = size(A, 1);
        D2 = D(1:sz, 1:sz);

        A = D2' * A * D;
        B = D2' * B * D;

        % Pack the rotated matricies into a single ABBA object
        if nargout == 2
          A = [ A B ; B A ];
        end
      elseif nargout ~= 4 && nargout ~= 1
        error('Insufficient number of output arguments');
      end
    end

    function [beam, D] = rotate(beam, varargin)
      %ROTATE apply the rotation matrix R to the beam coefficients
      %
      % [beam, D] = ROTATE(R) generates the wigner rotation matrix, D,
      % to rotate the beam.  R is the Cartesian rotation matrix
      % describing the rotation.  Returns the rotated beam and D.
      %
      % beam = ROTATE('wigner', D) applies the precomputed rotation.
      % This will only work if Nmax ~= 1.  If D is a cell array of
      % wigner matrices, generates a array of rotated beams.
      %
      % ROTATE(..., 'Nmax', nmax) specifies the Nmax for the
      % rotation matrix.  The beam Nmax is unchanged.  If nmax is smaller
      % than the beam Nmax, this argument is ignored.

      p = inputParser;
      p.addOptional('R', []);
      p.addParameter('Nmax', beam.Nmax);
      p.addParameter('wigner', []);
      p.parse(varargin{:});

      if ~isempty(p.Results.R) && isempty(p.Results.wigner)

        R = p.Results.R;

        % If no rotation, don't calculate wigner rotation matrix
        if sum(sum((eye(3) - R).^2)) < 1e-6
          D = eye(ott.utils.combined_index(p.Results.Nmax, p.Results.Nmax));
          return;
        end

        D = ott.utils.wigner_rotation_matrix(...
            max(beam.Nmax, p.Results.Nmax), R);
        beam = beam.rotate('wigner', D);

      elseif ~isempty(p.Results.wigner) && isempty(p.Results.R)
        
        if iscell(p.Results.wigner)
          
          ibeam = beam;
          beam = ott.Bsc();

          for ii = 1:numel(p.Results.wigner)
            sz = size(ibeam.a, 1);
            D2 = p.Results.wigner{ii}(1:sz, 1:sz);
            beam = beam.append(D2 * ibeam);
          end
        
        else

          sz = size(beam.a, 1);
          D2 = p.Results.wigner(1:sz, 1:sz);
          beam = D2 * beam;
          
        end

      else
        error('One of wigner or R must be specified');
      end
    end

    function [beam, D] = rotateX(beam, angle, varargin)
      %ROTATEX rotates the beam about the x-axis an angle in radians
      import ott.utils.*;
      [beam, D] = beam.rotate(rotx(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateY(beam, angle, varargin)
      %ROTATEX rotates the beam about the y-axis an angle in radians
      import ott.utils.*;
      [beam, D] = beam.rotate(roty(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateZ(beam, angle, varargin)
      %ROTATEX rotates the beam about the z-axis an angle in radians
      import ott.utils.*;
      [beam, D] = beam.rotate(rotz(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateXy(beam, anglex, angley, varargin)
      %ROTATEX rotates the beam about the x then y axes
      import ott.utils.*;
      [beam, D] = beam.rotate(roty(angley*180/pi)*rotx(anglex*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateXz(beam, anglex, anglez, varargin)
      %ROTATEX rotates the beam about the x then z axes
      import ott.utils.*;
      [beam, D] = beam.rotate(rotz(anglez*180/pi)*rotx(anglex*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateYz(beam, angley, anglez, varargin)
      %ROTATEX rotates the beam about the y then z axes
      import ott.utils.*;
      [beam, D] = beam.rotate(rotz(anglez*180/pi)*roty(angley*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateXyz(beam, anglex, angley, anglez, varargin)
      % Rotate the beam about the x, y then z axes
      %
      % [beam, D] = rorateXyz(anglex, angley, anglez, ...) additional
      % arguments are passed to beam.rotate.  Angles in radians.

      import ott.utils.*;
      [beam, D] = beam.rotate(rotz(anglez*180/pi)* ...
          roty(angley*180/pi)*rotx(anglex*180/pi), varargin{:});
    end

    function [a, b] = getCoefficients(beam, ci)
      %GETCOEFFICIENTS gets the beam coefficients
      %
      % ab = beam.getCoefficients() gets the beam coefficients packed
      % into a single vector, suitable for multiplying by a T-matrix.
      %
      % [a, b] = beam.getCoefficients() get the coefficients in two
      % beam vectors.
      %
      % beam.getCoefficients(ci) behaves as above but only returns
      % the requested beam cofficients a(ci) and b(ci).

      % If ci omitted, return all a and b
      if nargin == 1
        ci = 1:size(beam.a, 1);
      end

      a = beam.a(ci, :);
      b = beam.b(ci, :);

      if nargout == 1
        a = [a; b];
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

    function [moment, int, data] = intensityMoment(beam, varargin)
      % intensityMoment Calculate moment of beam intensity in the far-field
      %
      % [moment, int, data] = intensityMoment(...) integrates over the
      % incoming/outgoing field in the far-field to calculate the
      % moment of the intensity.  Also calculates the intensity.
      %
      % Can be used to calculate the force by comparing the outgoing component
      % of the incident beam with the total scattered beam.
      %
      % Optional named arguments:
      %   thetaRange   [min, max]  Range of angles from the pole to
      %      integrate over.  Default is 0 to pi (exclusive).
      %   saveData   bool  save data for repeated calculation (default: false)
      %   data       data  data saved for repeated calculation.
      %   ntheta     num   number of samples over 0 to pi theta range.
      %   nphi       num   number of samples over 0 to 2*pi phi range.

      % Parse inputs
      p = inputParser;
      p.addParameter('thetaRange', [0, pi]);
      p.addParameter('saveData', false);
      p.addParameter('ntheta', 100);
      p.addParameter('nphi', 100);
      p.addParameter('data', []);
      p.parse(varargin{:});

      % Regular beams have trivial solution
      if strcmpi(beam.basis, 'regular')
        moment = [0;0;0];
        int = 0;
        data = p.Results.data;
        warning('Regular wavefunctions go to zero in far-field');
        return;
      end

      % Setup the angular grid
      [theta, phi] = ott.utils.angulargrid(p.Results.ntheta, p.Results.nphi);
      dtheta = theta(2) - theta(1);
      dphi = phi(p.Results.ntheta+1) - phi(1);

      % Truncate the theta range
      keep = theta > p.Results.thetaRange(1) & theta < p.Results.thetaRange(2);
      theta = theta(keep);
      phi = phi(keep);

      uxyz = ott.utils.rtp2xyz([ones(size(theta)), theta, phi]).';
      
      % So integrals match sign convention used in ott.forcetorque
      uxyz(3, :) = -uxyz(3, :);

      % Calculate E-field in far-field
      [E, ~, data] = beam.farfield(theta, phi, ...
          'saveData', p.Results.saveData, ...
          'data', p.Results.data);

      % Calculate the irradiance
      Eirr = sum(abs(E).^2, 1);
      int = sum(Eirr .* sin(theta.') .* dtheta .* dphi, 2);

      % Calculate moment in Cartesian coordinates
      Eirr_xyz = uxyz .* Eirr;
      moment = sum(Eirr_xyz .* sin(theta.') .* dtheta .* dphi, 2);
    end

    function varargout = forcetorque(ibeam, other, varargin)
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

    function varargout = force(ibeam, other, varargin)
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
      
      % Combine incoherent beams
      if incN > 1
        fx = sum(reshape(fx, incN, []), 1);
        fy = sum(reshape(fy, incN, []), 1);
        fz = sum(reshape(fz, incN, []), 1);
      end
      
      % Ensure things are full
      fx = full(fx);
      fy = full(fy);
      fz = full(fz);

      % Package output
      if nargout == 3
        varargout{1:3} = {fx, fy, fz};
      else
        varargout{1} = [fx; fy; fz];
      end
    end

    function varargout = torque(ibeam, other, varargin)
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
        varargout{1} = [tx; ty; tz];
      end
    end

    function varargout = spin(ibeam, other, varargin)
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
        varargout{1} = [sx; sy; sz];
      end
    end

    function [sbeam, tbeam] = scatter(beam, tmatrix, varargin)
      % Calculate the beam scattered by a T-matrix
      %
      % Usage
      %   [sbeam, tbeam] = beam.scatter(tmatrix) scatters the beam
      %   returning the scattered beam ``sbeam`` and the unscattered
      %   but possibly translated beam ``tbeam`` truncated to
      %   ``tmatrix.Nmax + 1``.
      %
      % Optional named arguments
      %   - position (3xN numeric) translation applied to the beam
      %     before the beam is scattered by the particle.  Default: ``[]``.
      %
      %   - rotation (3x3N numeric) rotation spplied to the beam,
      %     calculates the scattered beam and applies the inverse rotation,
      %     effectively rotating the particle.  Default: ``[]``.
      %
      %   - combine (enum|empty) -- Beam combination method.  Can be
      %     'incoherent' (ignored), 'coherent', or empty (ignored).
      %
      % If both position and rotation are present, the translation is
      % applied first, followed by the rotation.
      % If both position and rotation are arrays, they must have the same
      % number of locations.

      dummy_type = 'total';
      tbeam = ott.optics.vswf.bsc.Scattered(beam, dummy_type);
      [sbeam, tbeam] = tbeam.calculateScatteredBeam(tmatrix, varargin{:});
    end

    function beam = mtimes(a,b)
      %MTIMES (op) divide the beam coefficients by a scalar
      %
      % Supports:
      %    - Scalar multiplication
      %    - Matrix multiplication of a and b vectors: A*a, a*A
      %    - Matrix multiplication of [a;b] vector: T*[a;b]

      if isa(a, 'ott.Bsc')
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

    function beam = plus(beam1, beam2)
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
    
    function beam = clearDz(beam)
      % Clear dz
      %
      % Useful when generating beams using translations.
      beam.dz = 0.0;
    end
  end

  methods (Hidden)
    function p = getBeamPower(beam)
      % get.power calculate the power of the beam
      p = full(sum(abs(beam.a).^2 + abs(beam.b).^2));
    end
    function beam = setBeamPower(beam, p)
      % set.power set the beam power
      beam = sqrt(p / beam.power) * beam;
    end

    function E = efieldInternal(beam, xyz)
      % Method used by efield(xyz)
      [E, ~, ~] = beam.emFieldXyz(xyz);
    end
    function H = hfieldInternal(beam, xyz)
      % Method used by hfield(xyz)
      [~, H, ~] = beam.emFieldXyz(xyz);
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
      sbeam = sbeam.totalField(ibeam);

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
    
    function [ibeam, sbeam, incN] = forcetorqueParser(ibeam, other, varargin)
      % Input parser for forcetorque and related methods
      
      ip = inputParser;
      ip.addParameter('position', []);  % Only used for T-matrix
      ip.addParameter('rotation', []);  % Only used for T-matrix
      ip.addParameter('combine', []);
      ip.parse(varargin{:});
      
      % Ensure we have a T-matrix
      if isa(other, 'ott.optics.vswf.tmatrix.Tmatrix')
        
        % Get the number of beams (for combination)
        % If coherent, combine before scattering
        Nbeams = ibeam.Nbeams;
        if strcmpi(ip.Results.combine, 'coherent')
          ibeam = sum(ibeam);
          Nbeams = 1;
        end
        
        [sbeam, ibeam] = ibeam.scatter(other, ...
            'position', ip.Results.position, ...
            'rotation', ip.Results.rotation);
      else
        sbeam = other;
        
        % Get the number of beams (for combination)
        assert(ibeam.Nbeams == 1 || ibeam.Nbeams == sbeam.Nbeams, ...
            'Number of incident and scattered beams must be 1 or matching');
        Nbeams = max(ibeam.Nbeams, sbeam.Nbeams);
      end

      % Handle the combine argument
      incN = 1;
      if isempty(ip.Results.combine) || Nbeams == 1
        % Nothing to do
        
      elseif strcmpi(ip.Results.combine, 'coherent')
        % Combine the beams before force calculation
        ibeam = sum(ibeam);
        sbeam = sum(sbeam);
        
      elseif strcmpi(ip.Results.combine, 'incoherent')
        % Return how many terms need to be combined later
        assert(numel(ibeam.Nbeams) == 1 ...
            || numel(ibeam.Nbeams) == numel(sbeam.Nbeams), ...
            'Number of incident and scattered beams must be 1 or match');
        incN = Nbeams;
        
      end
    end
  end
end
