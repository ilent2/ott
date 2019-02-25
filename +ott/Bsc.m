classdef Bsc
%Bsc abstract class representing beam shape coefficients
%
% Any units can be used for the properties as long as they are
% consistent in all specified properties.  Calculated quantities
% will have these units.
%
% Properties
%   a               Beam shape coefficients a vector
%   b               Beam shape coefficients b vector
%   type            Beam type (incident, scattered, total)
%   basis           VSWF beam basis (incoming, outgoing or regular)
%   Nmax            Truncation number for VSWF coefficients
%   power           Power of the beam [M*L^2/S^2]
%   Nbeams          Number of beams in this Bsc object
%   wavelength      Wavelength of beam [L]
%   speed           Speed of beam in medium [L/T]
%   omega           Angular frequency of beam [2*pi/T]
%   k_medium        Wavenumber in medium [2*pi/L]
%   dz              Absolute cumulative distance the beam has moved
%
% Methods
%   translateZ      Translates the beam along the z axis
%   translateXyz    Translation to xyz using rotations and z translations
%   translateRtp    Translation to rtp using rotations and z translations
%   farfield        Calculate fields in farfield
%   emFieldXyz      Calculate field values in cartesian coordinates
%   emFieldRtp      Calculate field values in spherical coordinates
%   getCoefficients Get the beam coefficients [a, b]
%   getModeIndices  Get the mode indices [n, m]
%   totalField      Calculate the total field reprsentation of the beam
%   scatteredField  Calcualte the scattered field representation of the beam
%   visualise       Generate a visualisation of the beam near-field
%   visualiseFarfield Generate a visualisation of the beam far-field
%   visualiseFarfieldSlice   Generate scattering slice at specific angle
%   visualiseFarfieldSphere  Generate spherical surface visualisation
%   intensityMoment Calculate moment of beam intensity in the far-field
%
% Static methods:
%   make_beam_vector    Convert output of bsc_* functions to beam coefficients
%
% Abstract methods:
%
% See also Bsc, ott.BscPmGauss, ott.BscPlane.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    a           % Beam shape coefficients a vector
    b           % Beam shape coefficients b vector

    omega       % Angular frequency of beam
    k_medium    % Wavenumber in medium

    dz          % Absolute cumulative distance the beam has moved

    % These can't be tracked using Matrix translation/rotations
    %offset      % Offset applied to beam using translate functions
    %direction   % Direction of beam applied using rotation functions
  end

  properties
    basis       % VSWF beam basis (incoming, outgoing or regular)
    type        % Beam type (incident, scattered, total)
  end

  properties (Dependent)
    Nmax        % Truncation number for VSWF coefficients
    power       % Power of the beam
    Nbeams      % Number of beams in this Bsc object

    wavelength  % Wavelength of beam
    speed       % Speed of beam in medium
  end

  methods (Abstract)
  end

  methods (Static)
    function [a, b, n, m] = make_beam_vector(a, b, n, m, Nmax)
      %MAKE_BEAM_VECTOR converts output of bsc_* functions to sparse vectors

      if isempty(n)
        error('OTT:BSC:make_beam_vector:no_modes', 'No modes');
      end

      if nargin < 5
        Nmax = max(n);
      end

      total_orders = ott.utils.combined_index(Nmax, Nmax);
      ci = ott.utils.combined_index(n, m);
      nbeams = size(a, 2);

      [ci, cinbeams] = meshgrid(ci, 1:nbeams);

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
    
    function data = GetVisualisationData(field_type, xyz, rtp, vxyz, vrtp)
      % Helper to generate the visualisation data output
      %
      % GetVisualisationData(field_type, xyz, rtp, vxyz, vrtp)
      % Takes a field_type string, the coordinates (either xyz or rtp),
      % and the data values (either vxyz or vrtp).
      %
      % vectors should all be Nx3 arrays to match inputs of
      % ott.utils.xyz2rtp and similar functions.
      %
      % field_type can be one of: 'irradiance', 'Et', 'E*',
      % 'Re(E*)', 'Abs(E*)' where the wildcard (*) must be one of
      % 'r', 't', 'p', 'x', 'y', 'z'.
      
      assert(size(xyz, 2) == 3 || size(xyz, 2) == 0, ...
        'xyz must be Nx3 matrix');
      assert(size(vxyz, 2) == 3 || size(vxyz, 2) == 0, ...
        'vxyz must be Nx3 matrix');
      assert(size(rtp, 2) == 3 || size(rtp, 2) == 0, ...
        'rtp must be Nx3 matrix');
      assert(size(vrtp, 2) == 3 || size(vrtp, 2) == 0, ...
        'vrtp must be Nx3 matrix');
      
      % Get the coordinates
      if isempty(xyz) && ~isempty(rtp)
        xyz = ott.utils.rtp2xyz(rtp);
      elseif isempty(rtp) && ~isempty(xyz)
        rtp = ott.utils.xyz2rtp(xyz);
      elseif isempty(rpt) && isempty(xyz)
        error('OTT:BSC:GetVisualisationData:no_coords', ...
          'Must supply coordinates');
      end
      
      % Get the data
      if isempty(vxyz) && ~isempty(vrtp)
        vxyz = ott.utils.rtpv2xyzv(vrtp, rtp);
      elseif isempty(vrtp) && ~isempty(vxyz)
        vrtp = ott.utils.xyzv2rtpv(vxyz, xyz);
      elseif isempty(vrtp) && isempty(vxyz)
        error('OTT:BSC:GetVisualisationData:no_data', ...
          'Must supply data');
      else
        error('OTT:BSC:GetVisualisationData:too_much_data', ...
          'Must supply only one data variable');
      end
      
      % Generate the requested field
      if strcmpi(field_type, 'irradiance')
        data = sqrt(sum(abs(vxyz).^2, 2));
        
      elseif strcmpi(field_type, 'Re(Er)')
        data = real(vrtp(:, 1));
      elseif strcmpi(field_type, 'Re(Et)')
        data = real(vrtp(:, 2));
      elseif strcmpi(field_type, 'Re(Ep)')
        data = real(vrtp(:, 3));
        
      elseif strcmpi(field_type, 'Re(Ex)')
        data = real(vxyz(:, 1));
      elseif strcmpi(field_type, 'Re(Ey)')
        data = real(vxyz(:, 2));
      elseif strcmpi(field_type, 'Re(Ez)')
        data = real(vxyz(:, 3));
        
      elseif strcmpi(field_type, 'Abs(Er)')
        data = abs(vrtp(:, 1));
      elseif strcmpi(field_type, 'Abs(Et)')
        data = abs(vrtp(:, 2));
      elseif strcmpi(field_type, 'Abs(Ep)')
        data = abs(vrtp(:, 3));
        
      elseif strcmpi(field_type, 'Abs(Ex)')
        data = abs(vxyz(:, 1));
      elseif strcmpi(field_type, 'Abs(Ey)')
        data = abs(vxyz(:, 2));
      elseif strcmpi(field_type, 'Abs(Ez)')
        data = abs(vxyz(:, 3));
        
      elseif strcmpi(field_type, 'Er')
        data = vrtp(:, 1);
      elseif strcmpi(field_type, 'Et')
        data = vrtp(:, 2);
      elseif strcmpi(field_type, 'Ep')
        data = vrtp(:, 3);
        
      elseif strcmpi(field_type, 'Ex')
        data = vxyz(:, 1);
      elseif strcmpi(field_type, 'Ey')
        data = vxyz(:, 2);
      elseif strcmpi(field_type, 'Ez')
        data = vxyz(:, 3);

      else
        error('OTT:BSC:GetVisualisationData:unknown_field_type', ...
          'Unknown field type value');
      end
      
    end
  end

  methods
    function beam = Bsc(a, b, basis, type, varargin)
      %BSC construct a new beam object
      %
      % beam = Bsc(a, b, basis, type, ...) constructs a new beam vector.
      % Useful if you have a specific set of a/b coefficients that you
      % want to wrap in a beam object.
      %
      % Optional named arguments:
      %    k_medium  n  Wavenumber in medium (default: 2*pi)
      %    omega     n  Angular frequency (default: 2*pi)
      %    dz        n  Initial displacement of the beam (default: 0)
      %    like    beam Construct this beam to be like another beam
      
      p = inputParser;
      p.addParameter('like', []);
      p.addParameter('k_medium', 2.0*pi);
      p.addParameter('omega', 2*pi);
      p.addParameter('dz', 0.0);
      p.parse(varargin{:});
      
      beam.dz = p.Results.dz;
      beam.k_medium = p.Results.k_medium;
      beam.omega = p.Results.omega;
      
      if ~isempty(p.Results.like)
        beam.omega = p.Results.like.omega;
        beam.k_medium = p.Results.like.k_medium;
        beam.dz = p.Results.like.dz;
      end

      if nargin ~= 0
        beam.a = a;
        beam.b = b;
        beam.basis = basis;
        beam.type = type;
      end
    end

    function beam = append(beam, other)
      % APPEND joins two beam objects together

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

      % SI-ify units of H
      H = H * -1i;
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
      p.addParameter('coord', 'cartesian');
      p.parse(varargin{:});

      % Scale the locations by the wave number (unitless coordinates)
      rtp(1, :) = rtp(1, :) * abs(beam.k_medium);

      % Get the indices required for the calculation
      [n,m]=ott.utils.combined_index(find(abs(beam.a)|abs(beam.b)));
      nm = [ n; m ];

      ci = ott.utils.combined_index(n, m);
      [a, b] = beam.getCoefficients(ci);

      % Calculate the fields
      [E, H, data] = ott.utils.emField(rtp.', beam.basis, nm, [a; b], ...
          'saveData', p.Results.saveData, ...
          'data', p.Results.data, ...
          'calcE', p.Results.calcE, 'calcH', p.Results.calcH);

      % Convert from spherical to Cartesian coordinates
      switch p.Results.coord
        case 'cartesian'
          E = ott.utils.rtpv2xyzv(E,rtp.');
          E(isnan(E)) = 0;
          H = ott.utils.rtpv2xyzv(H,rtp.');
          H(isnan(H)) = 0;

        case 'spherical'
          % Nothing to do

        otherwise
          error('Unknown coordinate system for output');
      end
      
      % Make the matrices 3xN
      E = E.';
      H = H.';
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

      rtp = ott.utils.xyz2rtp(xyz.');
      [E, H, data] = beam.emFieldRtp(rtp.', varargin{:});
    end
    
    function visualiseFarfieldSlice(beam, phi, varargin)
      % Generate a 2-D scattering plot of the far-field
      %
      % beam.visualiseFarfieldSlice(phi)
      
      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('normalise', false);
      p.addParameter('ntheta', 100);
      p.parse(varargin{:});
      
      ptheta = linspace(0, 2*pi, p.Results.ntheta);
      
      % TODO: Other field types

      [E, ~] = beam.farfield(ptheta, phi);
      I = sum(abs(E).^2, 1);
      
      if p.Results.normalise
        I = I ./ max(abs(I(:)));
      end

      polarplot(ptheta, I);
      
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
      dataout = beam.GetVisualisationData(p.Results.field, [], ...
        [itheta, iphi, ones(size(iphi))], [], ioutputE.');

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
        varargout{1} = im;
      elseif nargout == 2
        varargout{1} = im;
        varargout{2} = data;
      end
    end

    function varargout = visualise(beam, varargin)
      % Create a visualisation of the beam
      %
      % visualise(...) displays an image of the beam in the current
      % figure window.
      %
      % im = visualise(...) returns a image of the beam.
      %
      % Optional named arguments:
      %     'size'    [ x, y ]    Width and height of image
      %     'field'   type        Type of field to calculate
      %     'axis'    ax          Axis to visualise ('x', 'y', 'z') or
      %       a cell array with 2 or 3 unit vectors for x, y, [z].
      %     'offset'  offset      Plane offset along axis (default: 0.0)
      %     'range'   [ x, y ]    Range of points to visualise.
      %         Can either be a cell array { x, y }, two scalars for
      %         range [-x, x], [-y, y] or 4 scalars [ x0, x1, y0, y1 ].
      %     'mask'    func(xyz)   Mask function for regions to keep in vis

      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('size', [ 80, 80 ]);
      p.addParameter('axis', 'z');
      p.addParameter('offset', 0.0);
      p.addParameter('range', ...
          [1,1]*ott.utils.nmax2ka(beam.Nmax)/abs(beam.k_medium));
      p.addParameter('saveData', nargout == 2);
      p.addParameter('data', []);
      p.addParameter('mask', []);
      p.addParameter('showVisualisation', nargout == 0);
      p.parse(varargin{:});

      if iscell(p.Results.range)
        xrange = p.Results.range{1};
        yrange = p.Results.range{2};
        sz = [length(yrange), length(xrange)];
      elseif length(p.Results.range) == 2
        xrange = linspace(-1, 1, p.Results.size(1))*p.Results.range(1);
        yrange = linspace(-1, 1, p.Results.size(2))*p.Results.range(2);
        sz = p.Results.size;
      elseif length(p.Results.range) == 4
        xrange = linspace(-p.Results.range(1), p.Results.range(2), p.Results.size(1));
        yrange = linspace(-p.Results.range(3), p.Results.range(4), p.Results.size(2));
        sz = p.Results.size;
      else
        error('ott:Bsc:visualise:range_error', 'Incorrect number of range arguments');
      end
      [xx, yy, zz] = meshgrid(xrange, yrange, p.Results.offset);

      % Generate the xyz grid for the used requested plane
      if ischar(p.Results.axis)
        switch p.Results.axis
          case 'x'
            xyz = [zz(:), yy(:), xx(:)];
            alabels = {'Z', 'Y'};
          case 'y'
            xyz = [yy(:), zz(:), xx(:)];
            alabels = {'Z', 'X'};
          case 'z'
            xyz = [xx(:), yy(:), zz(:)];
            alabels = {'X', 'Y'};
          otherwise
            error('Unknown axis name specified');
        end
      elseif iscell(p.Results.axis)
        dir1 = p.Results.axis{1}(:);
        dir2 = p.Results.axis{2}(:);
        if numel(p.Results.axis) == 3
          dir3 = p.Results.axis{3}(:);
        else
          dir3 = cross(dir1(:), dir2(:));
        end
        
        alabels = {'Direction 1', 'Direction 2'};
        
        xyz = dir1.*xx(:).' + dir2.*yy(:).' + dir3.*zz(:).';
        xyz = xyz.';
      else
        error('axis must be character or cell array');
      end

      % Calculate the electric field
      [E, ~, data] = beam.emFieldXyz(xyz.', ...
          'saveData', p.Results.saveData, 'data', p.Results.data, ...
          'calcE', true', 'calcH', false);
        
      % Generate the requested field
      dataout = beam.GetVisualisationData(p.Results.field, ...
        xyz, [], E.', []);
      
      % Reshape the output
      imout = reshape(dataout, sz(2), sz(1));

      % Display the visualisation
      if p.Results.showVisualisation
        
        % Apply the mask
        if ~isempty(p.Results.mask)
          imout(p.Results.mask(xyz.')) = NaN;
        end
        
        imagesc(xrange, yrange, imout, 'AlphaData', ~isnan(imout));
        xlabel(alabels{1});
        ylabel(alabels{2});
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

    function p = get.power(beam)
      % get.power calculate the power of the beam
      p = full(sum(abs(beam.a).^2 + abs(beam.b).^2));
    end

    function beam = set.power(beam, p)
      % set.power set the beam power
      beam = sqrt(p / beam.power) * beam;
    end

    function lambda = get.wavelength(beam)
      % Get the beam wavelength
      lambda = 2*pi/beam.k_medium;
    end

    function beam = set.wavelength(beam, lambda)
      % Set the beam wavelength
      beam.k_medium = 2*pi/lambda;
    end

    function speed = get.speed(beam)
      % Get the speed of the beam in medium
      speed = beam.omega / beam.k_medium;
    end

    function beam = set.basis(beam, basis)
      % Set the beam type, checking it is a valid type first
      if ~any(strcmpi(basis, {'incoming', 'outgoing', 'regular'}))
        error('OTT:Bsc:set_basis:invalid_value', 'Invalid beam basis');
      end
      beam.basis = basis;
    end

    function beam = set.type(beam, type)
      % Set the beam type, checking it is a valid type first
      if ~any(strcmpi(type, {'incident', 'scattered', 'total', 'internal'}))
        error('OTT:Bsc:set_type:invalid_value', 'Invalid beam type');
      end
      beam.type = type;
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
      %TRANSLATEZ translate a beam along the z-axis
      %
      % TRANSLATEZ(z) translates by a distance z along the z axis.
      %
      % [beam, A, B] = TRANSLATEZ(z) returns the translation matrices
      % and the translated beam.  See also Bsc.TRANSLATE.
      %
      % [beam, AB] = TRANSLATEZ(z) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam.
      %
      % TRANSLATEZ(..., 'Nmax', Nmax) specifies the output beam Nmax.
      % Takes advantage of not needing to calculate a full translation matrix.

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
        if beam.dz > ott.utils.nmax2ka(beam.Nmax)/beam.k_medium
          warning('ott:Bsc:translateZ:outside_nmax', ...
              'Repeated translation of beam outside Nmax region');
        end
        beam.dz = beam.dz + abs(z);

        % Convert to beam units
        z = z * beam.k_medium / 2 / pi;

        ibeam = beam;
        beam = ott.Bsc();

        for ii = 1:numel(z)
          [A, B] = ott.utils.translate_z([p.Results.Nmax, ibeam.Nmax], z(ii));
          beam = beam.append(ibeam.translate(A, B));
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
      %TRANSLATEXYZ translate the beam given Cartesian coordinates
      %
      % TRANSLATEXYZ(xyz) translate the beam to locations given by
      % the xyz coordinates, where xyz is a 3xN matrix of coordinates.
      %
      % TRANSLATEXYZ(Az, Bz, D) translate the beam using
      % z-translation and rotation matricies.
      %
      % [beam, Az, Bz, D] = TRANSLATEXYZ(...) returns the z-translation
      % matrices, the rotation matrix D, and the translated beam.
      %
      % [beam, A, B] = TRANSLATEXYZ(...) returns the translation matrices
      % and the translated beam.
      %
      % [beam, AB] = TRANSLATEXYZ(...) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam.
      %
      % TRANSLATEXYZ(..., 'Nmax', Nmax) specifies the output beam Nmax.
      % Takes advantage of not needing to calculate a full translation matrix.

      p = inputParser;
      p.addOptional('opt1', []);    % xyz or Az
      p.addOptional('opt2', []);    % [] or Bz
      p.addOptional('opt3', []);    % [] or D
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      if ~isempty(p.Results.opt1) && isempty(p.Results.opt2) ...
          && isempty(p.Results.opt3)
        xyz = p.Results.opt1;
        rtp = ott.utils.xyz2rtp(xyz.').';
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
        beam = ott.Bsc();

        for ii = 1:numel(r)
          [tbeam, D] = ibeam.rotateYz(theta(ii), phi(ii), ...
              'Nmax', max(oNmax, ibeam.Nmax));
          [tbeam, A, B] = tbeam.translateZ(r(ii), 'Nmax', oNmax);
          beam = beam.append(tbeam.rotate('wigner', D'));
        end
      else
        dnmax = max(oNmax, beam.Nmax);
        D = eye(ott.utils.combined_index(dnmax, dnmax));

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
      % This will only work if Nmax ~= 1.
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

        sz = size(beam.a, 1);
        D2 = p.Results.wigner(1:sz, 1:sz);
        beam = D2 * beam;

      else
        error('One of wigner or R must be specified');
      end
    end

    function [beam, D] = rotateX(beam, angle, varargin)
      %ROTATEX rotates the beam about the x-axis an angle in radians
      [beam, D] = beam.rotate(rotx(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateY(beam, angle, varargin)
      %ROTATEX rotates the beam about the y-axis an angle in radians
      [beam, D] = beam.rotate(roty(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateZ(beam, angle, varargin)
      %ROTATEX rotates the beam about the z-axis an angle in radians
      [beam, D] = beam.rotate(rotz(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateXy(beam, anglex, angley, varargin)
      %ROTATEX rotates the beam about the x then y axes
      [beam, D] = beam.rotate(roty(angley*180/pi)*rotx(anglex*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateXz(beam, anglex, anglez, varargin)
      %ROTATEX rotates the beam about the x then z axes
      [beam, D] = beam.rotate(rotz(anglez*180/pi)*rotx(anglex*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateYz(beam, angley, anglez, varargin)
      %ROTATEX rotates the beam about the y then z axes
      [beam, D] = beam.rotate(rotz(anglez*180/pi)*roty(angley*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateXyz(beam, anglex, angley, anglez, varargin)
      % Rotate the beam about the x, y then z axes
      %
      % [beam, D] = rorateXyz(anglex, angley, anglez, ...) additional
      % arguments are passed to beam.rotate.  Angles in radians.

      [beam, D] = beam.rotate(rotz(anglez*180/pi)* ...
          roty(angley*180/pi)*rotx(anglex*180/pi), varargin{:});
    end

    function beam = totalField(beam, ibeam)
      % Calculate the total field representation of the beam
      %
      % total_beam = beam.totalField(incident_beam)
      
      switch beam.type
        case 'total'
          % Nothing to do
          
        case 'scattered'
          beam = 2*beam + ibeam;
          beam.type = 'total';
          
        case 'internal'
          error('Cannot convert from internal to total field');
          
        case 'incident'
          error('Cannot convert from incident to total field');
          
        otherwise
          error('Unknown beam type');
      end
    end

    function beam = scatteredField(beam, ibeam)
      % Calculate the scattered field representation of the beam
      %
      % scattered_beam = beam.totalField(incident_beam)
      
      switch beam.type
        case 'total'
          beam = 0.5*(beam - ibeam);
          beam.type = 'scattered';
          
        case 'scattered'
          % Nothing to do
          
        case 'internal'
          error('Cannot convert from internal to scattered field');
          
        case 'incident'
          error('Cannot convert from incident to total field');
          
        otherwise
          error('Unknown beam type');
      end
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

    function [sbeam, beam] = scatter(beam, tmatrix, varargin)
      %SCATTER scatter a beam using a T-matrix
      %
      % [sbeam, beam] = SCATTER(beam, tmatrix) scatters the beam
      % returning the scattered beam, sbeam, and the unscattered
      % but possibly translated beam truncated to tmatrix.Nmax + 1.
      %
      % SCATTER(..., 'position', xyz) applies a translation to the beam
      % before the beam is scattered by the particle.
      %
      % SCATTER(..., 'rotation', R) applies a rotation to the beam,
      % calculates the scattered beam and applies the inverse rotation,
      % effectively rotating the particle.

      % TODO: Support for multiple rotations or translations

      p = inputParser;
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.parse(varargin{:});

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
        maxNmax2 = min(maxNmax2, beam.Nmax);
      end

      % Ensure all T-matrices are the same size
      for ii = 1:numel(tmatrix)
        tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
      end

      % Apply translation to the beam
      if ~isempty(p.Results.position)

        % Requires scattered beam, convert if needed
        if ~strcmpi(tmatrix(1).type, 'scattered')
          maxNmax2 = min(maxNmax2, beam.Nmax);
          for ii = 1:numel(tmatrix)
            tmatrix(ii).type = 'scattered';
            tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
          end
        end

        % Apply translation
        % We need Nmax+1 terms for the force calculation
        beam = beam.translateXyz(p.Results.position, 'Nmax', maxNmax2+1);
      end

      % Apply rotation to the beam
      rbeam = beam;
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
      sbeam = ott.Bsc();
      for ii = 1:numel(tmatrix)
        sbeam = sbeam.append(tmatrix(ii).data * rbeam);
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
          sbeam.type = 'total';
          sbeam.basis = 'outgoing';
          
        case 'scattered'
          sbeam.type = 'scattered';
          sbeam.basis = 'outgoing';
          
        case 'internal'
          sbeam.type = 'internal';
          sbeam.basis = 'regular';
        
          % Wavelength has changed, update it
          sbeam.k_medium = tmatrix(1).k_particle;
          
        otherwise
          error('Unrecognized T-matrix type');
      end
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
    
    function beam = clearDz(beam)
      % Clear dz
      %
      % Useful when generating beams using translations.
      beam.dz = 0.0;
    end
  end
end
