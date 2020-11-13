function im = gaussfilt(im, filtSize, allwaysInternal)
% Applies Gaussian filtering to the input image
%
% If the image processing toolbox is installed, simply calls ``imgaussfilt``,
% otherwise uses ``conv2`` to apply a similar 2-D Gaussian convolution.
% Results wont be exactly the same.
%
% Usage
%   im = gaussfilt(im, filtSize)
%
%   im = gaussfilt(im, filtSize, alwaysInteranl)
%   Specify true as a third argument to always use the internal method.
%
% Parameters
%   - filtSize (numeric) -- Filter size in pixels.

% Copyright 2018-2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin == 2
  allwaysInternal = false;
end

if exist('imgaussfilt') == 2 && ~allwaysInternal
  im = imgaussfilt(im, filtSize);
else

  rng = ceil(2*filtSize);
  [X, Y] = meshgrid(-rng:rng, -rng:rng);
  G = exp(-(X.^2 + Y.^2)./(2*filtSize^2));
  im = conv2(im, G, 'same');

end

