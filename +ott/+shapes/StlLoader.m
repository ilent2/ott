classdef StlLoader < ott.shapes.TriangularMesh
% StlLoader load a shape from a STL file
%
% Properties:
%   filename   name of the file this object loaded
%   verts      (TriangularMesh) 3xN matrix of vertex locations
%   faces      (TriangularMesh) 3xN matrix of vertex indices describing faces
%   maxRadius  (Shape) maximum distance from shape origin
%   volume     (Shape) volume of shape
%
% Inherited methods:
%   writeWavefrontObj(shape, ...) write shape to Wavefront OBJ file.
%   insideXyz(shape, ...) determine if Cartesian point is inside shape.
%   voxels(shape, ...) xyz coordinates for voxels inside the shape.
%   surf(shape, ...) shape surface representation.
%
% See also StlLoader, ott.shapes.TriangularMesh, ott.shapes.WavefrontObj.
%
% This class uses a 3rd party STL file reader:
%  https://au.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader
% See tplicenses/stl_EricJohnson.txt for information about licensing.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    filename        % Name of the file this object loaded
  end

  methods (Hidden, Static)
    function [F, V, N] = stlbinary(M)
      % Function to read STL file
      %
      % This is from STL file reader version 1.2.0.0 (1.6 MB) by Eric Johnson
      % See tplicenses/stl_EricJohnson.txt for information about licensing.

      F = [];
      V = [];
      N = [];

      if length(M) < 84
        error('MATLAB:stlread:incorrectFormat', ...
              'Incomplete header information in binary STL file.');
      end

      % Bytes 81-84 are an unsigned 32-bit integer specifying the
      % number of faces that follow.
      numFaces = typecast(M(81:84),'uint32');
      %numFaces = double(numFaces);
      if numFaces == 0
        warning('MATLAB:stlread:nodata','No data in STL file.');
        return
      end

      T = M(85:end);
      F = NaN(numFaces,3);
      V = NaN(3*numFaces,3);
      N = NaN(numFaces,3);

      numRead = 0;
      while numRead < numFaces
        % Each facet is 50 bytes
        %  - Three single precision values specifying the face normal vector
        %  - Three single precision values specifying the first vertex (XYZ)
        %  - Three single precision values specifying the second vertex (XYZ)
        %  - Three single precision values specifying the third vertex (XYZ)
        %  - Two unused bytes
        i1    = 50 * numRead + 1;
        i2    = i1 + 50 - 1;
        facet = T(i1:i2)';

        n  = typecast(facet(1:12),'single');
        v1 = typecast(facet(13:24),'single');
        v2 = typecast(facet(25:36),'single');
        v3 = typecast(facet(37:48),'single');

        n = double(n);
        v = double([v1; v2; v3]);

				% Figure out where to fit these new vertices, and the face, in the
        % larger F and V collections.
        fInd  = numRead + 1;
        vInd1 = 3 * (fInd - 1) + 1;
        vInd2 = vInd1 + 3 - 1;

        V(vInd1:vInd2,:) = v;
        F(fInd,:)        = vInd1:vInd2;
        N(fInd,:)        = n;

        numRead = numRead + 1;
			end

    end
  end

  methods
    function shape = StlLoader(filename)
      % Construct a new shape from a STL file
      %
      % StlLoader(filename) loads the face and vertex information
      % contained in the file.
      %
      % Only supports binary STL files.
      %
      % This function uses 3rd party code,
      % see tplicenses/stl_EricJohnson.txt for licensing information.

      assert(nargin == 1, 'Filename not supplied');

      %
      % Begin third-party code
      %

      if ~exist(filename,'file')
        error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
               ', be sure to specify the full path to the file.'], filename);
      end

      fid = fopen(filename,'r');
      if ~isempty(ferror(fid))
          error(lasterror); %#ok
      end

      M = fread(fid,inf,'uint8=>uint8');
      fclose(fid);

      [faces,verts,~] = ott.shapes.StlLoader.stlbinary(M);

      %
      % End third-party code
      %

      % Convert to a mesh with unique vertices

      % Start by finding unique vertices
      [uverts,~,new_idx] = unique(verts,'rows');

      % Replace face indices with new indices
      faces = new_idx(faces);

      % Setup object
      shape = shape@ott.shapes.TriangularMesh(uverts.', faces.');
      shape.filename = filename;

    end
  end
end
