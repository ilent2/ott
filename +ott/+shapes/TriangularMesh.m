classdef TriangularMesh < ott.shapes.Shape
% TriangularMesh base class for triangular mesh objects (such as file loaders)
%
% Properties (read-only):
%   verts       3xN matrix of vertex locations
%   faces       3xN matrix of vertex indices describing faces
%
% Faces vertices should be ordered so normals face outwards for
% volume and inside functions to work correctly.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    verts           % Matrix of vertices in the object
    faces           % Matrix of faces in the object
  end

  methods
    function shape = TriangularMesh(verts, faces)
      % Construct a new triangular mesh representation
      %
      % TriangularMesh(verts, faces)
      %   verts       3xN matrix of vertex locations
      %   faces       3xN matrix of vertex indices describing faces
      %
      % Faces vertices should be ordered so normals face outwards for
      % volume and inside functions to work correctly.

      shape = shape@ott.shapes.Shape();

      % Verify size of inputs
      assert(size(verts, 1) == 3, 'Verts must be matrix of 3xN');
      assert(size(faces, 1) == 3, 'Faces must be matrix of 3xN');

      shape.verts = verts;
      shape.faces = faces;

      % Verify we have sufficient verts for faces
      if max(shape.faces(:)) > numel(shape.verts)
        error('faces matrix refers to non-existent vertices');
      end

    end

    function r = get_maxRadius(shape)
      % Calculate the maximum distance from the shape origin
      r = max(sum(shape.verts.^2, 1).^0.5);
    end

    function totalVolume = get_volume(shape)
      % Calculate the volume of the shape
      %
      % This function is based on a surface triangulation volume code
      % version 1.0.0.0 (1.43 KB) by Krishnan Suresh
      % matlabcentral/fileexchange/26982-volume-of-a-surface-triangulation
      % See tplicenses/stl_KrishnanSuresh.txt for information about licensing.

      p = shape.verts;
      t = shape.faces;

      % Compute the vectors d13 and d12
      d13= [(p(1,t(2,:))-p(1,t(3,:))); (p(2,t(2,:))-p(2,t(3,:))); ...
          (p(3,t(2,:))-p(3,t(3,:)))];
      d12= [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); ...
          (p(3,t(1,:))-p(3,t(2,:)))];

      cr = cross(d13,d12,1);%cross-product (vectorized)

      area = 0.5*sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);% Area of each triangle
      totalArea = sum(area);

      crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
      zMean = (p(3,t(1,:))+p(3,t(2,:))+p(3,t(3,:)))/3;
      nz = -cr(3,:)./crNorm;% z component of normal for each triangle

      volume = area.*zMean.*nz; % contribution of each triangle
      totalVolume = sum(volume);%divergence theorem

    end

    function b = inside(shape, radius, theta, phi)
      % Determine if spherical point point is inside shape
      %
      % b = inside(shape, radius, theta, phi) determine if the
      % point described by radius, theta (polar), phi (azimuthal)
      % is inside the shape.

      [x, y, z] = ott.utils.rtp2xyz(radius(:), theta(:), phi(:));
      b = shape.insideXyz(x, y, z);
    end

    function b = insideXyz(shape, x, varargin)
      % INSIDEXYZ determine if Cartesian point is inside the shape
      %
      % b = inside(shape, x, y, z) determine if the Cartesian point
      % [x, y, z] is inside the star shaped object.
      %
      % b = insideXyz(shape, xyz) as above, but using a 3xN matrix of
      % [x; y; z] positions.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.
      %
      % See also INSIDE.

      p = inputParser;
      p.addOptional('y', [], @isnumeric);
      p.addOptional('z', [], @isnumeric);
      p.addParameter('origin', 'world');
      p.parse(varargin{:});
      
      if isempty(p.Results.y) && isempty(p.Results.z)
        y = x(2, :);
        z = x(3, :);
        x = x(1, :);
      elseif ~isempty(p.Results.y) && ~isempty(p.Results.z)
        x = x(:);
        y = p.Results.y(:);
        z = p.Results.z(:);
        [x, y, z] = ott.utils.matchsize(x, y, z);
      else
        error('Must suply either 3xN matrix or x, y and z');
      end

      % Translate to shape origin
      if strcmpi(p.Results.origin, 'world')
        x = x - shape.position(1);
        y = y - shape.position(2);
        z = z - shape.position(3);
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''');
      end

      % Using a third-party function for insidexyz
      b = ott.utils.inpolyhedron(shape.faces.', shape.verts.', ...
          [x(:),y(:),z(:)]);
    end

    function surf(shape, varargin)
      % SURF generate a visualisation of the shape
      %
      % SURF(...) displays a visualisation of the shape in the current figure.
      %
      % Optional named arguments:
      %   offset   [x;y;z]   offset for location of surface
      %   rotation   mat     rotation matrix to apply to surface
      %   surfoptions   {varargin} options to be passed to surf.

      p = inputParser;
      p.addParameter('surfoptions', {});
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.addParameter('axes', []);
      p.parse(varargin{:});
      
      X = shape.verts(1, :);
      Y = shape.verts(2, :);
      Z = shape.verts(3, :);
      
      % Apply a rotation to the surface
      if ~isempty(p.Results.rotation)
        XYZ = [X; Y; Z];
        XYZ = p.Results.rotation * XYZ;
        X = XYZ(1, :);
        Y = XYZ(2, :);
        Z = XYZ(3, :);
      end

      % Apply a translation to the surface
      if ~isempty(p.Results.position)
        X = X + p.Results.position(1);
        Y = Y + p.Results.position(2);
        Z = Z + p.Results.position(3);
      end
      
      if nargout == 0 || ~isempty(p.Results.axes)
        
        % Place the surface in the specified axes
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = axes();
        end
        
        % Create the surface in the background, copy it to the axes
        f = figure('visible','off');
        h = trisurf(shape.faces.', X, Y, Z);
        c = copyobj(h, our_axes);
        set(c, p.Results.surfoptions{:});
        c.Visible = 'on';
        close(f);
      end

    end

    function writeWavefrontObj(shape, filename)
      % Write internal representation of shape to Wavefront OBJ file
      %
      % writeWavefrontObj(filename) writes the shape to the given file.
      shape.writeWavefrontObj_helper(filename, shape.verts, shape.faces);
    end
  end
end
