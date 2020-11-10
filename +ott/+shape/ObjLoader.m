classdef ObjLoader < ott.shape.TriangularMesh
% Load a shape from a Wavefront OBJ file.
% Inherits from :class:`TriangularMesh`.
%
% Properties
%   - filename      -- Filename for loaded OBJ file
%
% The file format is described at
% `<https://en.wikipedia.org/wiki/Wavefront_.obj_file>`_

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    filename        % Name of the file this object loaded
  end

  methods
    function shape = ObjLoader(filename)
      % Construct a new shape from a Wavefront OBJ file
      %
      % Loads the face and vertex information
      % contained in the file.  Faces are converted to triangles.
      %
      % Usage
      %   shape = ObjLoader(filename)

      assert(nargin == 1, 'Filename not supplied');

      fp = fopen(filename, 'r');

      verts = zeros(3, 0);
      faces = zeros(3, 0);

      % Read lines of file
      fline = fgets(fp);
      while fline ~= -1

        switch fline(1:2)
          case 'v '

            % Vert lines should have 3 or more coordinates: x y z [w]
            v = sscanf(fline(3:end), '%f');

            if length(v) < 3
              warning('ott:shape:ObjLoader:vertex_with_two_points', ...
                  'Ignoring vertex with insufficient coords');
            else
              verts(:, end+1) = v(1:3); %#ok<AGROW>
            end

          case 'f '

            % Read face indices, they can have multiple formats:
            %    f 1 2 3
            %    f 1/2 2/2 3/2
            %    f 1//2 2//2 3//2
            f1 = sscanf(fline(3:end), '%f');
            f2 = sscanf(fline(3:end), '%f%*[^ ]');

            % Determine which pattern matched
            if length(f1) > length(f2)
              f = f1;
            else
              f = f2;
            end

            if length(f) < 3
              warning('ott:shape:ObjLoader:face_with_too_few_verts', ...
                  'Ignoring face with fewer than 3 verts');
            else
              for ii = 3:length(f)
                faces(:, end+1) = [f(1); f(ii-1); f(ii)]; %#ok<AGROW>
              end
            end

          otherwise
            % Nothing to do, ignore line
        end

        % Get next line of file
        fline = fgets(fp);
      end

      % Close file
      fclose(fp);

      % Setup object
      shape = shape@ott.shape.TriangularMesh(verts, faces);
      shape.filename = filename;

    end
  end
end
