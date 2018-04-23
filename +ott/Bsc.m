classdef Bsc
%Bsc abstract class representing beam shape coefficients
%
% Bsc properties:
%   a               Beam shape coefficients a vector
%   b               Beam shape coefficients b vector
%
% Bsc methods:
%   translateZ      Translates the beam along the z axis
%   translateXyz    Translation to xyz using rotations and z translations
%   translateRtp    Translation to rtp using rotations and z translations
%   farfield        Calculate fields in farfield
%   emfieldXyz      Calculate fields at specified locations
%   set.Nmax        Resize the beam shape coefficient vectors
%   get.Nmax        Get the current size of the beam shape coefficient vectors
%   getCoefficients Get the beam coefficients [a, b]
%   getModeIndices  Get the mode indices [n, m]
%   power           Calculate the power of the beam
%
% Static methods:
%   make_beam_vector    Convert output of bsc_* functions to beam coefficients
%
% Abstract methods:
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    a           % Beam shape coefficients a vector
    b           % Beam shape coefficients b vector
    type        % Coefficient type (incomming, outgoing or regular)

    k_medium    % Wavenumber in medium

    dz          % Absolute cumulative distance the beam has moved

    % These can't be tracked using Matrix translation/rotations
    %offset      % Offset applied to beam using translate functions
    %direction   % Direction of beam applied using rotation functions
  end

  properties (Dependent)
    Nmax        % Size of beam vectors
    power       % Power of the beam
    beams       % Number of beams in this Bsc object
  end

  methods (Abstract)
  end

  methods (Static)
    function [a, b, n, m] = make_beam_vector(a, b, n, m, Nmax)
      %MAKE_BEAM_VECTOR converts output of bsc_* functions to sparse vectors

      if isempty(n)
        error('No modes');
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
  end

  methods
    function beam = Bsc(a, b, type, varargin)
      %BSC construct a new beam object

      if nargin ~= 0
        beam.a = a;
        beam.b = b;
        beam.type = type;
      end

      beam.dz = 0.0;
    end

    function beam = append(beam, other)
      % APPEND joins two beam objects together

      beam.Nmax = max(beam.Nmax, other.Nmax);
      other.Nmax = beam.Nmax;
      beam.a = [beam.a, other.a];
      beam.b = [beam.b, other.b];
    end

    function [E, H] = farfield(beam, theta, phi)
      %FARFIELD finds far field at locations thieta, phi.

      [theta,phi] = ott.utils.matchsize(theta,phi);

      [theta_new,~,indY]=unique(theta);
      [phi_new,~,indP]=unique(phi);

      Etheta=zeros(length(theta),1);
      Ephi=zeros(length(theta),1);

      Htheta=zeros(length(theta),1);
      Hphi=zeros(length(theta),1);

      if strcmp(beam.type, 'incomming')

        a = beam.a;
        b = beam.b;
        p = zeros(size(beam.a));
        q = zeros(size(beam.b));

      elseif strcmp(beam.type, 'outgoing')

        a = zeros(size(beam.a));
        b = zeros(size(beam.a));
        p = beam.a;
        q = beam.b;

      else

        % TODO: Can we convert from regular to incomming + outgoing?
        error('Unsupported beam type');

      end

      a = ott.utils.threewide(a);
      b = ott.utils.threewide(b);
      p = ott.utils.threewide(p);
      q = ott.utils.threewide(q);

      [n,m]=ott.utils.combined_index(find(abs(beam.a)|abs(beam.b)));

      for nn = 1:max(n)

        vv=find(n==nn);
        [Y,Ytheta,Yphi] = ott.utils.spharm(nn,m(vv), ...
            theta_new,zeros(size(theta_new)));

        [M,PHI]=meshgrid(m(vv),phi_new);

        expimphi=repmat(exp(1i*M.*PHI),[1,3]);

        %this makes the vectors go down in m for n.
        % has no effect if old version code.
        Nn = 1/sqrt(nn*(nn+1));

        for ii=1:length(vv)
          index=nn*(nn+1)+m(vv(ii));

          TEMP=Nn*(((1i)^(nn+1) * a(index)+(-1i)^(nn+1) ...
              * p(index)).*Yphi(:,ii)+((1i)^nn * b(index) ...
              +(-1i)^nn * q(index)).*Ytheta(:,ii));
          Etheta=Etheta+TEMP(indY).*expimphi(indP,ii);

          TEMP=Nn*(-((1i)^(nn+1) * a(index)+(-1i)^(nn+1) ...
              * p(index)).*Ytheta(:,ii)+((1i)^nn * b(index) ...
              +(-1i)^nn * q(index)).*Yphi(:,ii));
          Ephi=Ephi+TEMP(indY).*expimphi(indP,ii);

          if nargout>1
            TEMP=Nn*(((1i)^(nn+1) * b(index)+(-1i)^(nn+1) ...
                * q(index)).*Yphi(:,ii)+((1i)^nn * a(index) ...
                +(-1i)^nn * p(index)).*Ytheta(:,ii));
            Htheta=Htheta+TEMP(indY).*expimphi(indP,ii);

            TEMP=Nn*(-((1i)^(nn+1) * b(index)+(-1i)^(nn+1) ...
                * q(index)).*Ytheta(:,ii)+((1i)^nn * a(index) ...
                +(-1i)^nn * p(index)).*Yphi(:,ii));
            Hphi=Hphi+TEMP(indY).*expimphi(indP,ii);
          end
        end
      end

      E=[zeros(size(Etheta)),Etheta,Ephi].';
      H=[zeros(size(Htheta)),Htheta,Hphi].';

      % SI-ify units of H
      H = H * -1i;
    end

    function [E, H] = emFieldXyz(beam, xyz)
      %EMFIELDXYZ calculates the E and H field at specified locations

      kxyz = xyz * beam.k_medium;

      [n,m]=ott.utils.combined_index(find(abs(beam.a)|abs(beam.b)));
      nm = [ n; m ];

      if strcmp(beam.type, 'incomming')
        S = ott.electromagnetic_field_xyz(kxyz.', nm, beam, [], []);
        E = S.Eincident.';
        H = S.Hincident.';
      elseif strcmp(beam.type, 'outgoing')
        S = ott.electromagnetic_field_xyz(kxyz.', nm, [], beam, []);
        E = S.Escattered.';
        H = S.Hscattered.';
      elseif strcmp(beam.type, 'regular')
        S = ott.electromagnetic_field_xyz(kxyz.', nm, [], [], beam);
        E = S.Einternal.';
        H = S.Hinternal.';
      else
        error('Invalid beam type');
      end
    end

    function p = get.power(beam)
      % get.power calculate the power of the beam
      p = sum(abs(beam.a).^2 + abs(beam.b).^2);
    end

    function beam = set.power(beam, p)
      % set.power set the beam power
      beam = sqrt(p / beam.power) * beam;
    end

    function beam = set.type(beam, type)
      % Set the beam type, checking it is a valid type first
      if ~any(strcmpi(type, {'incomming', 'outgoing', 'regular'}))
        error('Invalid beam type');
      end
      beam.type = type;
    end

    function beams = get.beams(beam)
      % get.beams get the number of beams in this object
      beams = size(beam.a, 2);
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

        beam.a = beam.a(1:total_orders);
        beam.b = beam.b(1:total_orders);

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
        beam.a = sparse(arow_index,acol_index,aa,total_orders,1);
        beam.b = sparse(brow_index,bcol_index,ba,total_orders,1);
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

        [A, B] = ott.utils.translate_z([p.Results.Nmax, beam.Nmax], z);
      else
        error('Wrong number of arguments');
      end

      % Apply the translation
      beam = beam.translate(A, B);

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
        rtp = ott.utils.xyz2rtp(xyz(:).');
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

      % Handle input arguments
      if ~isempty(p.Results.opt1) && isempty(p.Results.opt2) ...
          && isempty(p.Results.opt3)

        % Assume first argument is rtp coordinates
        r = p.Results.opt1(1);
        theta = p.Results.opt1(2);
        phi = p.Results.opt1(3);

      elseif ~isempty(p.Results.opt1) && ~isempty(p.Results.opt2) ...
          && ~isempty(p.Results.opt3)

        % Rotation/translation is already computed, apply it
        A = p.Results.opt1;
        B = p.Results.opt2;
        D = p.Results.opt3;
        beam = D * beam;
        beam = beam.translate(A, B);

        % The beam might change size, so readjust D to match
        sz = size(A, 1);
        D2 = D(1:sz, 1:sz);
        beam = D2' * beam;
        return;
      else
        error('Not enough input arguments');
      end

      % Only do the rotation if we need it
      if theta ~= 0 || phi ~= 0
        [beam, D] = beam.rotateYz(theta, phi);
        [beam, A, B] = beam.translateZ(r, 'Nmax', p.Results.Nmax);

        % The beam might change size, so readjust D to match
        sz = size(A, 1);
        D2 = D(1:sz, 1:sz);
        beam = D2' * beam;
      else
        D = eye(size(beam.a, 1));
        [beam, A, B] = beam.translateZ(r, 'Nmax', p.Results.Nmax);
      end

      % Rotate the translation matricies
      if nargout <= 3

        % The beam might change size, so readjust D to match
        sz = size(A, 1);
        D2 = D(1:sz, 1:sz);

        A = D2' * A * D;
        B = D2' * B * D;

        % Pack the rotated matricies into a single ABBA object
        if nargout == 2
          A = [ A B ; B A ];
        end
      end
    end

    function [beam, D] = rotate(beam, R)
      %ROTATE apply the rotation matrix R to the beam coefficients

      % If no rotation, don't calculate wigner rotation matrix
      if sum(sum((eye(3) - R).^2)) < 1e-6
        D = eye(size(beam.a, 1));
        return;
      end

      D = ott.utils.wigner_rotation_matrix(beam.Nmax, R);
      beam = D * beam;
    end

    function [beam, D] = rotateX(beam, angle)
      %ROTATEX rotates the beam about the x-axis an angle in radians
      [beam, D] = beam.rotate(rotx(angle*180/pi));
    end

    function [beam, D] = rotateY(beam, angle)
      %ROTATEX rotates the beam about the y-axis an angle in radians
      [beam, D] = beam.rotate(roty(angle*180/pi));
    end

    function [beam, D] = rotateZ(beam, angle)
      %ROTATEX rotates the beam about the z-axis an angle in radians
      [beam, D] = beam.rotate(rotz(angle*180/pi));
    end

    function [beam, D] = rotateXy(beam, anglex, angley)
      %ROTATEX rotates the beam about the x then y axes
      [beam, D] = beam.rotate(roty(angley*180/pi)*rotx(anglex*180/pi));
    end

    function [beam, D] = rotateXz(beam, anglex, anglez)
      %ROTATEX rotates the beam about the x then z axes
      [beam, D] = beam.rotate(rotz(anglez*180/pi)*rotx(anglex*180/pi));
    end

    function [beam, D] = rotateYz(beam, angley, anglez)
      %ROTATEX rotates the beam about the y then z axes
      [beam, D] = beam.rotate(rotz(anglez*180/pi)*roty(angley*180/pi));
    end

    function [beam, D] = rotateXyz(beam, anglex, angley, anglez)
      %ROTATEX rotates the beam about the x, y then z axes
      [beam, D] = beam.rotate(rotz(anglez*180/pi)* ...
          roty(angley*180/pi)*rotx(anglex*180/pi));
    end

    function beam = outgoing(beam, ibeam)
      %TOOUTGOING calculate the outgoing beam
      if strcmp(beam.type, 'outgoing')
        % Nothing to do
      elseif strcmp(beam.type, 'regular')
        beam = 2*beam + ibeam;
        beam.type = 'outgoing';
      else
        error('Unable to convert incomming beam to outgoing beam');
      end
    end

    function beam = regular(beam, ibeam)
      %TOREGULAR calculate regular beam
      if strcmp(beam.type, 'outgoing')
        beam = 0.5*(beam - ibeam);
        beam.type = 'regular';
      elseif strcmp(beam.type, 'regular')
        % Nothing to do
      else
        error('Unable to convert incomming beam to outgoing beam');
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

    function [sbeam, beam] = scatter(beam, tmatrix, varargin)
      %SCATTER scatter a beam using a T-matrix
      %
      % [sbeam, beam] = SCATTER(beam, tmatrix) scatters the beam
      % returning the scattered beam, sbeam, and the unscattered
      % but possibly translated/rotated beam.
      %
      % SCATTER(..., 'position', xyz) applies a translation to the beam.
      % SCATTER(..., 'rotation', R) applies a rotation to the beam.

      p = inputParser;
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.parse(varargin{:});

      % Apply translation to the beam (requires scattered T-matrix)
      if ~isempty(p.Results.position)
        tmatrix.type = 'scattered';
        beam = beam.translateXyz(p.Results.position, 'Nmax', tmatrix.Nmax(2));
      end

      % Apply rotation to the beam
      if ~isempty(p.Results.rotation)
        beam = beam.rotate(p.Results.rotation);
      end

      % Ensure the Nmax for the inner dimension matches
      if strcmpi(tmatrix.type, 'scattered')
        newNmax = min(tmatrix.Nmax(2), beam.Nmax);
        tmatrix = tmatrix.set_Nmax([tmatrix.Nmax(1), newNmax], ...
            'powerloss', 'ignore');
        beam = beam.set_Nmax(newNmax, 'powerloss', 'ignore');
      else
        tmatrix = tmatrix.set_Nmax([tmatrix.Nmax(1), beam.Nmax], ...
            'powerloss', 'ignore');
      end

      % Calculate the resulting beam
      sbeam = tmatrix.data * beam;

      % Assign a type to the resulting beam
      if strcmp(tmatrix.type, 'total')
        sbeam.type = 'outgoing';
      elseif strcmp(tmatrix.type, 'scattered')
        sbeam.type = 'regular';
      elseif strcmp(tmatrix.type, 'internal')
        sbeam.type = 'regular';
      else
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
  end
end
