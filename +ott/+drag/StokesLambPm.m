classdef StokesLambPm < ott.drag.Stokes
% Calculate drag coefficients using Lamb series and point matching.
%
% Based on code by Lachlan Gibson

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function drag = simple(shape, varargin)
      % Calculate the drag using the Stokes sphere approximation.
      %
      % Usage:
      %   drag = ott.drag.Sphere.simple(shape, ...) construct a new
      %   drag tensor for the given ott.shapes.Shape object.
      %   Shape must implement a maxRadius method, the result is used
      %   as the sphere radius.
      %
      %   drag = ott.drag.Sphere.simple(name, parameters, ...) constructs
      %   a new shape described by name and parameters.
      %   See ott.shapes.Shape.simple for supported shapes.
      %   Constructed shape must have a maxRadius method.
      %
      % See ott.drag.Sphere/Sphere for optional arguments.

      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('parameters', []);
      p.addParameter('viscosity', 1.0);
      p.parse(varargin{:});

      % Get a shape object from the inputs
      if ischar(shape) && ~isempty(p.Results.parameters)
        shape = ott.shapes.Shape.simple(shape, p.Results.parameters);
        varargin = varargin(2:end);
      elseif ~isa(shape, 'ott.shapes.Shape') || ~isempty(p.Results.parameters)
        error('Must input either Shape object or string and parameters');
      end

      % Check the particle is star shaped
      if ~isa(shape, 'ott.shapes.StarShape')
        error('Only star shaped particles supported for now');
      end

      % Calculate drag
      drag = ott.drag.StokesLambPm(rtp, p.Results.viscosity);
    end
  end

  methods (Hidden)
    function forcetorque = calculate_column(rtp, v, eta, sorder, rmax)
        % Prepare coefficient values
        vals=[...
          Comp_Vals('Spherical',rtp.',v(:,1),1,eta,sorder,rmax),...
          Comp_Vals('Spherical',rtp.',v(:,2),2,eta,sorder,rmax),...
          Comp_Vals('Spherical',rtp.',v(:,3),3,eta,sorder,rmax)];

        [sol,NM,RMSE]=solve_sys(vals,[0 sorder],rmax);

        % Force
        force=4*pi*...
          [sol(NM(:,1)==1&NM(:,2)==1&NM(:,3)==1&NM(:,4)==1&NM(:,5)==2),...
          sol(NM(:,1)==1&NM(:,2)==1&NM(:,3)==1&NM(:,4)==2&NM(:,5)==2),...
          -sol(NM(:,1)==1&NM(:,2)==0&NM(:,3)==1&NM(:,4)==1&NM(:,5)==2)];

        % Torque
        torque=8*pi*eta*...
          [sol(NM(:,1)==1&NM(:,2)==1&NM(:,3)==3&NM(:,4)==1&NM(:,5)==2),...
          sol(NM(:,1)==1&NM(:,2)==1&NM(:,3)==3&NM(:,4)==2&NM(:,5)==2),...
          -sol(NM(:,1)==1&NM(:,2)==0&NM(:,3)==3&NM(:,4)==1&NM(:,5)==2)];

        forcetorque = [force(:); torque(:)];
    end
  end

  methods
    function obj = StokesLambPm(rtp, eta, sorder)

      % Calculate maximum radius
      rmax = max(rtp(1, :));

      xyz = ott.utils.rtp2xyz(rtp);

      forward = zeros(6, 6);

      % Build up 6 columns of drag tensor
      % TODO: Symmetries?
      vel = {[1,0,0], [0,1,0], [0,0,1]};
      for ii = 1:6

        switch
          case 1
            v = repmat([1,0,0],length(R),1);
          case 2
            v = repmat([0,1,0],length(R),1);
          case 3
            v = repmat([0,0,1],length(R),1);
          case 4
            v = [zeros(size(xyz(1, :))), -xyz(3, :), xyz(2, :)].';
          case 5
            v = [-xyz(3, :), zeros(size(xyz(1, :))), xyz(1, :)].';
          case 6
            v = [-xyz(2, :), xyz(1, :), zeros(size(xyz(1, :)))].';
        end

        forward(:, ii) = obj.calculate_column(rtp, v, eta, sorder, rmax);

      end

      obj = obj@Stokes('forward', forward);
    end
  end
end
