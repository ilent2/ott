classdef Shape
%Shape abstract class for optical tweezers toolbox shapes
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties
  end

  methods (Static)

    function shape = simple(name, parameters)
			%SIMPLE constructs a shape for a bunch of simple particles
      % This method replaces shapesurface.m from OTTv1.
      %
      % Supported shapes [parameters]:
      %   'sphere'          Sphere [ radius ]
      %   'cylinder'        z-axis aligned cylinder [ radius height ]
      %   'ellipsoid'       Ellipsoid [ a b c]
      %   'superellipsoid'  Superellipsoid [ a b c e n ]
      %   'cone-tipped-cylinder'      [ radius height cone_height ]
      %   'cube'            Cube [ width ]
      %   'axisym'          Axis-symetric particle [ rho(:) z(:) ]

      switch lower(name)
        case 'sphere'
          shape = ott.shapes.Sphere(parameters(:));
        case 'cylinder'
          shape = ott.shapes.Cylinder(parameters(:));
        case 'ellipsoid'
          shape = ott.shapes.Ellipsoid(parameters(:));
        case 'superellipsoid'
          shape = ott.shapes.Superellipsoid(parameters(:));
        case 'cone-tipped-cylinder'

          radius = parameters(1);
          height = parameters(2);
          cone_height = parameters(3);

          z = [ height/2 + cone_height, height/2, ...
              -height/2, -height/2 - cone_height];
          rho = [ 0.0, radius, radius, 0.0 ];

          shape = ott.shapes.AxisymLerp(rho, z);
        case 'cube'
          shape = ott.shapes.Cube(parameters(:));
        case 'axisym'
          shape = ott.shapes.AxisymLerp(parameters(:));
        otherwise
          error('Unsupported simple particle shape');
      end
    end
  end

  methods
  end
end
