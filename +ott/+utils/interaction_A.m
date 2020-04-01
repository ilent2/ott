function A = interaction_A(k,r,varargin)
% Calculate the interaction matrix
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
p.parse(varargin{:});

if ~isempty(p.Results.alpha) && ~isempty(p.Results.inv_alpha)
  error('Either alpha or inv_alpha must be supplied or none, not both');
end

assert(ismatrix(r) && size(r, 2) == 3, 'voxels must be a Nx3 matrix');

N = size(r, 1);

% Generate the matrix of off-diagonal elements
A = zeros(3*N,3*N);
for j=1:N
  A((1:3) + 3*(j-1),:) = ott.utils.calc_Aj(k,r,[],j,false);
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
  A = A + blkdiag(Ac{:});
end

