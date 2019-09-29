function [fx,fy,fz,tx,ty,tz,sx,sy,sz]=forcetorque(ibeam, sbeam, varargin)
% FORCETORQUE calculate force/torque/spin in a 3D orthogonal space
% If the beam shape coefficients are in the original coordinates,
% this outputs the force, torque and spin in 3D carteisan coordinates.
%
% Units are beam power.  Force results should be multipled by n/c
% and torque results multiplied by 1/omega, assuiming the beam coefficients
% already have the correct units for power.
%
% [fxyz,txyz,sxyz] = FORCETORQUE(ibeam, sbeam) calculates the force,
% torque and spin using the incident beam, ibeam, and the scattered
% beam, sbeam.
%
% Output is stored in [3, 1] column vectors.  If torque or spin are
% omitted, only force or force/torque are calculated.
%
% FORCETORQUE(ibeam, T, 'position', position) first applies a translation
% to the beam.  position can be a 3xN array, resulting in multiple
% force/torque calculations for each position.
%
% FORCETORQUE(ibeam, T, 'rotation', rotation) effectively applies a
% rotation to the particle by first applying the rotation to the beam,
% scattering the beam by the T-matrix and applying the inverse rotation
% to the scattered beam.  rotation can be a 3x3N array, resulting in
% multiple calculations.
%
% If both position and rotation are present,
% the translation is applied first, followed by the rotation.
% If both position and rotation are arrays, they must have the same
% number of locations (N) or a single location (N=1).
%
% ibeam can contain multiple beams.  If multiple beams are present,
% the outputs are [3, nlocations, nbeams] arrays unless the
% coherent argument is set to true, in which case the beams are added
% after translation.
%
% [fx,fy,fz,tx,ty,tz,sx,sy,sz] = FORCETORQUE(...) unpacks the
% force/torque/spin into separate output arguments.
%
% This uses mathematical result of Farsund et al., 1996, in the form of
% Chricton and Marsden, 2000, and our standard T-matrix notation S.T.
% E_{inc}=sum_{nm}(aM+bN);
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott.warning('internal');

fx=0;
fy=0;
fz=0;
tx=0;
ty=0;
tz=0;
sx=0;
sy=0;
sz=0;

% Parse optional inputs
p = inputParser;
p.addParameter('position', []);
p.addParameter('rotation', []);
p.addParameter('coherent', false);
p.addParameter('progress_callback', []);
p.parse(varargin{:});

% Check sizes of inputs
assert(size(p.Results.position, 1) == 3 ...
    || numel(p.Results.position) == 0, ...
    'position must either be a empty array or 3xN array');
assert(size(p.Results.rotation, 1) == 3 ...
    || numel(p.Results.rotation) == 0, ...
    'rotation must either be a empty array or 3x3N array');

if isa(sbeam, 'ott.Tmatrix')

  % Rename T-matrix
  T = sbeam;
  nparticles = numel(T);

  npositions = max(1, size(p.Results.position, 2));
  nrotations = max(1, size(p.Results.rotation, 2)/3);

  if npositions ~= 1 && nrotations ~= 1 && npositions ~= nrotations
    error('OTT:forcetorque:nlocations', ...
      'Number of positions/rotations should be equal or 1');
  end

  nlocations = max([npositions, nrotations]);

  % Ensure all T-matricies have appropriate type
  for ii = 1:nparticles
    T(ii).type = 'scattered';
  end
  
  nbeams = ibeam.Nbeams;
  if p.Results.coherent
    nbeams = 1;
  end

  % Preallocate output
  f = zeros(3*numel(T), nlocations, nbeams);
  t = f;
  s = f;

  for ii = 1:nlocations
    
    % Output the progress
    if ~isempty(p.Results.progress_callback)
      p.Results.progress_callback((ii-1)/nlocations);
    end

    position = [];
    if ~isempty(p.Results.position)
      if npositions == 1
        position = p.Results.position;
      else
        position = p.Results.position(:, ii);
      end
    end

    rotation = [];
    if ~isempty(p.Results.rotation)
      if nrotations == 1
        rotation = p.Results.rotation;
      else
        rotation = p.Results.rotation(:, 3*(ii-1) + (1:3));
      end
    end

    % Calculate the scattered beams and translated beam
    [sbeam, tbeam] = ibeam.scatter(T, ...
        'position', position, 'rotation', rotation);
      
    % If beams are coherent, combine them
    if p.Results.coherent
      sbeam = sbeam.mergeBeams();
      tbeam = tbeam.mergeBeams();
    end

    % Calculate force
    [fl,tl,sl] = ott.forcetorque(tbeam, sbeam);

    % Unpack the calculated force
    f(:, ii, :) = reshape(fl, 3*numel(T), 1, nbeams);
    t(:, ii, :) = reshape(tl, 3*numel(T), 1, nbeams);
    s(:, ii, :) = reshape(sl, 3*numel(T), 1, nbeams);
  end
  
  % Squeeze to 2-D array if Nbeams is 1
  f = squeeze(f);
  t = squeeze(t);
  s = squeeze(s);

  % Move output to appropriate locations
  if nargout > 3
    fx = f(1*(1:nparticles), :);
    fy = f(2*(1:nparticles), :);
    fz = f(3*(1:nparticles), :);
    tx = t(1*(1:nparticles), :);
    ty = t(2*(1:nparticles), :);
    tz = t(3*(1:nparticles), :);
    sx = s(1*(1:nparticles), :);
    sy = s(2*(1:nparticles), :);
    sz = s(3*(1:nparticles), :);
  else
    fx = f;
    fy = t;
    fz = s;
  end

  ott.warning('external');
  return;
end

% Check the number of beams in each input
if ibeam.Nbeams ~= sbeam.Nbeams && ibeam.Nbeams ~= 1 && sbeam.Nbeams ~= 1
  error('Beam objects must contain same number of beams or 1 beam');
end

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

Az=m./n./(n+1).*imag(-(a).*conj(b)+conj(q).*(p)); %this has the correct sign... farsund. modes match.
Bz=1./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ... %.*n
    .*imag(anp1.*conj(a)+bnp1.*conj(b)-(pnp1).*conj(p) ...
    -(qnp1).*conj(q)); %this has the correct sign... farsund. modes match.

fz=2*sum(Az+Bz);

Axy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)).*(conj(pmp1).*q - conj(amp1).*b - conj(qmp1).*p + conj(bmp1).*a); %this has the correct sign... farsund. modes match.
Bxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ... %sqrt(n.*)
    ( sqrt((n+m+1).*(n+m+2)) .* ( p.*conj(pnp1mp1) + q.* conj(qnp1mp1) -a.*conj(anp1mp1) -b.*conj(bnp1mp1)) + ... %this has the correct sign... farsund. modes match.
    sqrt((n-m+1).*(n-m+2)) .* (pnp1mm1.*conj(p) + qnp1mm1.*conj(q) - anp1mm1.*conj(a) - bnp1mm1.*conj(b)) ); %this has the correct sign... farsund. modes match.

fxy=sum(Axy+Bxy);
fx=real(fxy);
fy=imag(fxy);

if nargout > 1
    tz=sum(m.*(a.*conj(a)+b.*conj(b)-p.*conj(p)-q.*conj(q))); %this has the correct sign... farsund. modes match.
    
    txy=sum(sqrt((n-m).*(n+m+1)).*(a.*conj(amp1)+b.*conj(bmp1)-p.*conj(pmp1)-q.*conj(qmp1))); %this has the correct sign... farsund. modes match.
    tx=real(txy);
    ty=imag(txy);
    
    if nargout > 2
        Cz=m./n./(n+1).*(-(a).*conj(a)+conj(q).*(q)-(b).*conj(b)+conj(p).*(p));
        Dz=-2./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ...
            .*real(anp1.*conj(b)-bnp1.*conj(a)-(pnp1).*conj(q) ...
            +(qnp1).*conj(p));
        
        sz = sum(Cz+Dz);
        
        Cxy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)).*(conj(pmp1).*p - conj(amp1).*a + conj(qmp1).*q - conj(bmp1).*b);
        Dxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ...
            ( (sqrt((n+m+1).*(n+m+2)) .* ( p.*conj(qnp1mp1) - q.* conj(pnp1mp1) -a.*conj(bnp1mp1) +b.*conj(anp1mp1))) + ...
            (sqrt((n-m+1).*(n-m+2)) .* (pnp1mm1.*conj(q) - qnp1mm1.*conj(p) - anp1mm1.*conj(b) + bnp1mm1.*conj(a))) );
        
        sxy=sum(Cxy+Dxy);
        sy=real(sxy);
        sx=imag(sxy);
    end
end

if nargout <= 3
    fx=full([fx;fy;fz]);
    fy=full([tx;ty;tz]);
    fz=full([sx;sy;sz]);
end

ott.warning('external');
