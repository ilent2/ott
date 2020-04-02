function A = interaction_A(k,r,varargin)
% Calculate the interaction matrix
%
% TODO: Should this be part of TmatrixDda?
%
% A = interaction_A(k_medium, voxels, alpha) calculate the off-diagonal
% terms and use the polarisability, alpha, to form the diagonal.
% voxels must be a Nx3 array of Cartesian voxel coordinates.
%
% A = interaction_A(k_medium, voxels, 'inv_alpha', inv_alpha) as above,
% but uses the inverse polarisability.
%
% alpha and inv_alpha can be scalar and 3x1 vectors or 3x3 matrix for
% homogeneous isotropic and birefringent material.  For inhomogeneous
% materials, use N, 3xN (or 3N vector) or 3x3N matrices respectively.
% If alpha or inv_alpha are not supplied, the diagonal is left empty.
%
% Based on OMG/Code/DDA/interaction_A

% Copyright 2018 Isaac Lenton

p = inputParser;
p.addOptional('alpha', []);
p.addParameter('inv_alpha', []);
p.addParameter('z_rotational_symmetry', 1);
p.parse(varargin{:});

if ~isempty(p.Results.alpha) && ~isempty(p.Results.inv_alpha)
  error('Either alpha or inv_alpha must be supplied or none, not both');
end

assert(p.Results.z_rotational_symmetry >= 1, ...
  'Only discrete rotational symmetry supported for now');

assert(ismatrix(r) && size(r, 2) == 3, 'voxels must be a Nx3 matrix');

N = size(r, 1);

% Pre-compute cart2sph transformations for each dipole
cart2sph = [];
sph2cart = [];
if p.Results.z_rotational_symmetry > 1
  
  % Pre-allocate memory
  cart2sph = zeros(3, 3*N);
  sph2cart = cart2sph;
  
  % Compute
  for ii = 1:N
    [~, theta, phi] = ott.utils.xyz2rtp(r(ii, :));

    cart2sph(:, (1:3) + (ii-1)*3) = ott.TmatrixDda.cart2sph_mat(...
        theta, phi);
  end
end

% Generate the matrix of off-diagonal elements
A = zeros(3*N,3*N, p.Results.z_rotational_symmetry);
for m = 1:p.Results.z_rotational_symmetry
  
  % Pre-compute sph2cart transformations for each dipole
  if m > 1
    for ii = 1:N
      [~, theta, phi] = ott.utils.xyz2rtp(r(ii, :));

      sph2cart(:, (1:3) + (ii-1)*3) = ott.TmatrixDda.sph2cart_mat(...
          theta, phi + 2*pi*(m-1)/p.Results.z_rotational_symmetry);
    end
  end
  
  for j=1:N
    A((1:3) + 3*(j-1),:, m) = calc_Aj(k,r,j, p.Results.z_rotational_symmetry, m, ...
      cart2sph, sph2cart);
  end
end

% Put the inverse polarisability elements into usable form
inv_alpha = [];
if ~isempty(p.Results.alpha)
  val = p.Results.alpha;
  sz = size(val);
  if all(sz == [1, 1])
    inv_alpha = repmat([1./val, 0, 0; 0, 1./val, 0; 0, 0, 1./val], 1, N);
  elseif all(sz == [3, 3])
    inv_alpha = repmat(inv(val), 1, N);
  elseif all(sz == [3, 1])
    inv_alpha = repmat(diag(1.0./val), 1, N);
  elseif numel(val) == N
    inv_alpha = zeros(3, 3*N);
    inv_alpha(1, 1:3:end) = 1.0./val;
    inv_alpha(2, 2:3:end) = 1.0./val;
    inv_alpha(3, 3:3:end) = 1.0./val;
  elseif all(sz == [3, N])
    inv_alpha = zeros(3, 3*N);
    inv_alpha(1, 1:3:end) = 1.0./val(1:3:end);
    inv_alpha(2, 2:3:end) = 1.0./val(2:3:end);
    inv_alpha(3, 3:3:end) = 1.0./val(3:3:end);
  elseif all(sz == [3, 3*N])
    inv_alpha = zeros(3, 3*N);
    for ii = 1:N
      inv_alpha(:, (1:3) + 3*(ii-1)) = inv(val(:, (1:3) + 3*(ii-1)));
    end
  else
    error('Unsupported size of alpha');
  end
elseif ~isempty(p.Results.inv_alpha)
  val = p.Results.inv_alpha;
  sz = size(val);
  if all(sz == [1, 1])
    inv_alpha = repmat([val, 0, 0; 0, val, 0; 0, 0, val], 1, N);
  elseif all(sz == [3, 3])
    inv_alpha = repmat(val, 1, N);
  elseif all(sz == [3, 1])
    inv_alpha = repmat(diag(val), 1, N);
  elseif numel(val) == N
    inv_alpha = zeros(3, 3*N);
    inv_alpha(1, 1:3:end) = val;
    inv_alpha(2, 2:3:end) = val;
    inv_alpha(3, 3:3:end) = val;
  elseif all(sz == [3, N])
    inv_alpha = zeros(3, 3*N);
    inv_alpha(1, 1:3:end) = val(1:3:end);
    inv_alpha(2, 2:3:end) = val(2:3:end);
    inv_alpha(3, 3:3:end) = val(3:3:end);
  elseif all(sz == [3, 3*N])
    inv_alpha = val;
  else
    error('Unsupported size of inv_alpha');
  end
end

% Put the inverse polarisability on the diagonal
if ~isempty(inv_alpha)
  Ac = mat2cell(inv_alpha, 3, repmat(3,1,N));
  A(:, :, 1) = A(:, :, 1) + blkdiag(Ac{:});
end

end

function Aj = calc_Aj(k_medium, r, j, z_rotation, m, cart2sph, sph2cart)
% calculates a 3 X 3N block comprising N number of 3 X 3 Green's tensors

  % The following method is about twice as slow as the previous
  % implementation (for no rot. sym.), but it has a lot of similarity
  % with the F matrix and perhaps we can write a nice & optimal version
  % later!!!
  
  % Get voxel location for this quadrant
  if m == 1
    our_xyz = r;
  else
    rotphi = (m-1) * 2*pi/z_rotation;
    rotM = [cos(rotphi) -sin(rotphi) 0; sin(rotphi) cos(rotphi) 0; 0 0 1];
    our_xyz = (rotM * (r.')).';
  end

  n_dipoles = size(r, 1);

  r_vec = our_xyz - r(j, :);
  r_jk = vecnorm(r_vec, 2, 2);
  
  % Remove diagonal term, only applies when not calculating
  % elements with rotational symmetry.  (avoids nan)
  if m == 1
    r_jk(j) = 1;
  end
  
  r_hat = r_vec./r_jk;
  
  Aj = zeros(3, 3*n_dipoles);
  
  for ii = 1:n_dipoles
    
    % Skip self-interaction term
    if ii == j && m ~= 1
      continue;
    end
    
    % Calculate coordinate rotation
    if m == 1
      M = 1;
    else
      M = sph2cart(:, (1:3) + (ii-1)*3) * cart2sph(:, (1:3) + (ii-1)*3);
    end
    
    rr = r_hat(ii, :).'*r_hat(ii, :);
    
    Aj(:, (1:3) + 3*(ii-1)) = exp(1i*k_medium*r_jk(ii, :))/r_jk(ii, :)*...
      (k_medium^2*(rr - eye(3)) + (1i*k_medium*r_jk(ii, :)-1)/r_jk(ii, :)^2*(3*rr - eye(3))) * M;
    
  end
  
end
