classdef StokesLambNn < ott.drag.StokesStarShaped
% Calculate stokes drag for star shaped particles using pre-trained NN.
% Inherits from :class:`StokesStarShaped`.
%
% Uses the neural network from Lachlan Gibson's PhD thesis (2016).
%
% Properties
%   - inverse     -- Calculated from `forward`
%   - forward     -- Drag tensor calculated using point matching.
%   - viscosity   -- Viscosity of medium
%   - shape       -- A star shaped particle describing the geometry
%
% Static methods
%   - SphereGrid  -- Generate grid of angular points for NN
%   - LoadNetwork -- Load the neural network
%
% Additional methods/properties inherited from :class:`Stokes`.

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (Dependent)
    forwardInternal     % Calculate forward drag tensor
  end

  properties (Access=protected, Hidden)
    network
  end

  methods (Static)
    function [phi, theta] = SphereGrid()
      % Generate grid of angular points for evaluating neural network
      %
      % Usage
      %   [phi, theta] = StokesLambNn.SphereGrid()
      %
      % Generates 100 angular coordinates for evaluating the NN.
      % Based on Gibson's Sphere_Grid code.

      % Number of points
      n = 100;

      lin = linspace(-1,1,n); % Linear parameterisation
      theta = pi/2+asin(lin); % Polar angle
      a = sqrt(pi*n); % Length of Spiral is approximately 2a
      phi=theta * a; % Azimuthal angle
    end

    function net = LoadNetwork()
      % Load the neural network from a file
      %
      % Usage
      %   net = StokesLambNn.LoadNetwork()

      % Get current folder/package
      folder = fileparts(which('ott.drag.StokesLambNn'));

      % Load the network
      fdata = load(fullfile(folder, 'private', 'GibsonStar30Ffnet'), 'net');
      net = fdata.net;
    end
  end

  methods
    function drag = StokesLambNn(varargin)
      % Construct star-shaped drag method using Neural network approach.
      %
      % Usage
      %   drag = StokesLambNn(shape, ...)
      %
      % Parameters
      %   - shape (ott.shapes.Shape) -- Star shaped particle.
      %
      % Parameters can also be passed as named arguments.
      % Additional parameters passed to base.

      % Only need constructor for doc/help
      drag = drag@ott.drag.StokesStarShaped(varargin{:});

      % Load network on construction (rather than each time we need it)
      % in order to avoid thrashing the file system
      drag.network = drag.LoadNetwork();
    end
  end

  methods % Getters/setters
    function D = get.forwardInternal(drag)

      % Evaluate shape points
      [phi, theta] = drag.SphereGrid();
      radii = drag.shape.starRadii(theta, phi);

      % Scale by mean radius
      meanRadius = mean(radii);
      radii = radii ./ meanRadius;

      % Evaluate network
      T = drag.network(radii(:));

      % Unpack drag
      D = zeros(6, 6);
      for ii=1:6
        ind=22-(7-ii)*(8-ii)/2;
        D(ii,ii:end)=T(ind:ind+6-ii).';
        D(ii:end,ii)=T(ind:ind+6-ii);
      end

      % Scale by viscosity
      D = -D .* drag.viscosity;

      % Add factor for cross-terms
      D(1:3, 4:6) = D(1:3, 4:6) * (4/3);

      % Scale by radius
      D(1:3, 1:3) = D(1:3, 1:3) .* meanRadius;
      D(1:3, 4:6) = D(1:3, 4:6) .* meanRadius.^2;
      D(4:6, 1:3) = D(4:6, 1:3) .* meanRadius.^2;
      D(4:6, 4:6) = D(4:6, 4:6) .* meanRadius.^3;

    end
  end
end
