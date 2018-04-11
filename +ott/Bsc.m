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

    dz          % Distance the beam has been translated
  end

  properties (Dependent)
    Nmax        % Current size of beam vectors
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

      a = sparse(ci, 1, a, total_orders, 1);
      b = sparse(ci, 1, b, total_orders, 1);

      [n, m] = ott.utils.combined_index(1:Nmax^2+2*Nmax);
      n = n.';
      m = m.';
    end

    function k_medium = parser_k_medium(p)
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

    function nmax = get.Nmax(beam)
      %get.Nmax calculates Nmax from the current size of the beam coefficients
      nmax = ott.utils.combined_index(length(beam.a));
    end

    function beam = set.Nmax(beam, nmax)
      %set.Nmax resizes the beam vectors
      total_orders = ott.utils.combined_index(nmax, nmax);
      if length(beam.a) > total_orders

        amagA = full(sum(sum(abs(beam.a).^2)));
        bmagA = full(sum(sum(abs(beam.b).^2)));

        beam.a = beam.a(1:total_orders);
        beam.b = beam.b(1:total_orders);

        amagB = full(sum(sum(abs(beam.a).^2)));
        bmagB = full(sum(sum(abs(beam.b).^2)));

        aapparent_error = abs( amagA - amagB )/amagA;
        bapparent_error = abs( bmagA - bmagB )/bmagA;

        warning_error_level = 1e-6;
        if aapparent_error > warning_error_level || ...
            bapparent_error > warning_error_level
          warning('ott:Bsc:setNmax:truncation', ...
              ['Apparent errors of ' num2str(aapparent_error) ...
                  ', ' num2str(bapparent_error) ]);
        end
      elseif length(beam.a) < total_orders
        [arow_index,acol_index,aa] = find(beam.a);
        [brow_index,bcol_index,ba] = find(beam.b);
        beam.a = sparse(arow_index,acol_index,aa,total_orders,1);
        beam.b = sparse(brow_index,bcol_index,ba,total_orders,1);
      end
    end

    function beam = translateZ(beam, z)
      %TRANSLATEZ translate a beam along the z-axis

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
      
      [A, B] = ott.utils.translate_z(beam.Nmax, z);
      beam = [ A B ; B A ] * beam;
    end

    function beam = translateXyz(beam, x, y, z)
      %TRANSLATEXYZ translate the beam given Cartesian coordinates
      if nargin == 2
        rtp = ott.utils.xyz2rtp(x(:).');
      else
        rtp = ott.utils.xyz2rtp(x, y, z);
      end
      beam = beam.translateRtp(rtp);
    end

    function beam = translateRtp(beam, r, theta, phi)
      %TRANSLATERTP translate the beam given spherical coordinates
      if nargin == 2
        theta = r(2);
        phi = r(3);
        r = r(1);
      end

      % Only do the rotation if we need it
      if theta ~= 0 || phi ~= 0
        [beam, D] = beam.rotateYz(theta, phi);
        beam = beam.translateZ(r);
        beam = D' * beam;
      else
        beam = beam.translateZ(r);
      end
    end

    function [beam, D] = rotate(beam, R)
      %ROTATE apply the rotation matrix R to the beam coefficients
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

    function beam = toOutgoing(beam, ibeam)
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

    function beam = toRegular(beam, ibeam)
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
      [n, m] = ott.utils.combined_index([1:length(beam.a)].');
      if nargout == 1
        n = [n; m];
      end
    end

    function p = power(beam)
      %POWER calculate the power of the beam
      p = sqrt(sum(abs(beam.a).^2 + abs(beam.b).^2));
    end

    function beam = mrdivide(beam,o)
      %MRDIVIDE (op) divide the beam coefficients by a scalar
      beam.a = beam.a / o;
      beam.b = beam.b / o;
    end

    function beam = scatter(beam, tmatrix)
      %SCATTER scatter a beam using a T-matrix

      % Ensure the Nmax for the inner dimension matches
      tmatrix.Nmax = [tmatrix.Nmax, beam.Nmax];

      % Calculate the resulting beam
      beam = tmatrix.data * beam;

      % Assign a type to the resulting beam
      if strcmp(tmatrix.type, 'total')
        beam.type = 'regular';
      elseif strcmp(tmatrix.type, 'scattered')
        beam.type = 'outgoing';
      elseif strcmp(tmatrix.type, 'internal')
        beam.type = 'regular';
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
