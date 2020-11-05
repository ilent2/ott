function angle_rad = na2angle(na, refractive_index)
% Calculates the half-angle for a numerical aperture value.
%
% Calculates the angle using::
%
%   angle = asin(na ./ refractive_index)
%
% Usage
%   angle_rad = na2angle(na, refractive_index)
%
% Parameters
%   - na (numeric) -- Numerical aperture.
%
%   - refractive_index (numeric) -- Refractive index.  This is typically
%     the refractive index of the medium or objective (depending on the
%     application).  Default: ``1.0`` (i.e., vacuum or air).

if nargin == 1
  refractive_index = 1;
end

assert(isnumeric(refractive_index) && isscalar(refractive_index) ...
    && refractive_index > 0 && isfinite(refractive_index), ...
    'refractive index must be finite positive numeric scalar');
assert(isnumeric(na) && isscalar(na) && na >= 0 && na <= refractive_index, ...
    'na must be numeric scalar between 0 and refractive_index');

angle_rad = asin(na ./ refractive_index);

end

