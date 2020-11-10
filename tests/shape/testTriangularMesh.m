function tests = testTriangularMesh
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testConstructor(testCase)

  verts = randn(3, 3);
  faces = [1; 2; 3];

  shape = ott.shape.TriangularMesh(verts, faces);

  testCase.verifyEqual(shape.verts, verts, 'verts');
  testCase.verifyEqual(shape.faces, faces, 'faces');
  testCase.verifyEqual(shape.boundingBox, ...
    [min(verts, [], 2), max(verts, [], 2)], 'bb');
  
  % This should filter the bad face
  shape = ott.shape.TriangularMesh(verts, [faces, [1;1;3]]);
  testCase.verifyEqual(shape.verts, verts, 'bad face verts');
  testCase.verifyEqual(shape.faces, faces, 'bad face faces');
  
  % Scale
  shape = shape * 2;
  testCase.verifyEqual(shape.verts, verts.*2, 'scale');

end

function testSubdivide(testCase)

  shape = ott.shape.Cube();
  shape = ott.shape.TriangularMesh(shape);
  
  test = shape.subdivide();

end

function testIntersect(testCase)

  cube = ott.shape.Cube();
  shape = ott.shape.PatchMesh(cube);  % Uses TriangularMesh internally
  
  x0 = [0.1;0;-10];
  x1 = [0.0;0;10];
  
  testCase.verifyEqual(shape.intersect(x0, x1), ...
    cube.intersect(x0, x1), 'isect');
  
  shapeInts = shape.intersectAll(x0, x1, 'removeNan', true);
  cubeInts = cube.intersectAll(x0, x1, 'removeNan', true);
  testCase.verifyEqual(shapeInts, cubeInts, 'isectAll');
  

end

function testInsideXyz(testCase)

  cube = ott.shape.Cube();
  cube = ott.shape.TriangularMesh(cube);
  testCase.assertClass(cube, 'ott.shape.TriangularMesh', 'type');
  radius = cube.maxRadius;

  % Choose three points inside the shape and one outside
  btarget = [true, true, true, false].';
  x = [0.5*radius.*rand(1, 3), 4.0];
  y = [0.5*radius.*rand(1, 3), 4.0];
  z = [0.5*radius.*rand(1, 3), 4.0];

  xyz = [x(:), y(:), z(:)].';

  b = cube.insideXyz(xyz);
  testCase.verifyEqual(b, btarget, 'insideXyz');

end

function testStarRadiiCube(testCase)

  width = 1;
  cube = ott.shape.Cube(width);
  cube = ott.shape.TriangularMesh(cube);
  cube.starShaped = true;
  
  theta = [0, pi/2, pi, pi/2];
  phi = [0, pi/2, 0, pi/4];
  radii = cube.starRadii(theta, phi);
  
  testCase.verifyEqual(radii, [1,1,1, sqrt(2)]*width/2, ...
      'AbsTol', 1e-15, 'width');

end

function testStarRadiiSphere(testCase)

  radius = 2;
  sphere = ott.shape.Sphere(radius);
  sphere = ott.shape.TriangularMesh(sphere);
  sphere.starShaped = true;
  
  theta = randn(10, 10);
  phi = randn(10, 10);
  radii = sphere.starRadii(theta, phi);
  
  testCase.verifyEqual(radii, radius*ones(size(theta)), ...
      'AbsTol', 3.1e-2, 'radius');

end

function testSmallMesh(testCase)

  T =   [0    0.1653    0.3307    0.4960
         0    0.1653    0.3307    0.4960
         0    0.1653    0.3307    0.4960
         0    0.1653    0.3307    0.4960
         0    0.1653    0.3307    0.4960
         0    0.1653    0.3307    0.4960];
       
  P = [0         0         0         0
    0.3142    0.3142    0.3142    0.3142
    0.6283    0.6283    0.6283    0.6283
    0.9425    0.9425    0.9425    0.9425
    1.2566    1.2566    1.2566    1.2566
    1.5708    1.5708    1.5708    1.5708];
  
  R = [1.1931    0.8432    1.1010    0.9227
    1.1931    0.9945    0.8238    0.9222
    1.1931    1.0331    1.0363    0.8031
    1.1931    0.9401    0.9826    0.9420
    1.1931    0.9446    1.0438    1.0328
    1.1931    1.1624    1.0197    1.1617];
  
  [X, Y, Z] = ott.utils.rtp2xyz(R, T, P);
  
  shape = ott.shape.PatchMesh.FromSurfMatrix(X, Y, Z);
  shape.starShaped = true;
  
  radii = shape.starRadii(T, P);
  testCase.verifyEqual(sum(isnan(radii(:))), 0, 'all');
  
end

