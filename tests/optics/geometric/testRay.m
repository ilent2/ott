function tests = testRay
  % Test cases based on example_ray.m script from OTGO
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function test1dArray(testCase)

  num_rays = 5;
  origin = zeros(3, num_rays);
  direction = rand(3, num_rays);

  vectors = ott.utils.Vector(origin, direction);
  power = ones(1, num_rays);
  polarisation = cross(vectors, ott.utils.Vector(origin, ones(size(origin))));

  rays = ott.optics.geometric.Ray(vectors, power, polarisation);

  testCase.verifyEqual(size(rays), [num_rays, 1], ...
    'Size of rays incorrect');

end

function test2dArray(testCase)

  xrays = 5;
  yrays = 3;

  origin = zeros(3, yrays, xrays);
  direction = rand(3, yrays, xrays);

  vectors = ott.utils.Vector(origin, direction);
  power = ones(1, yrays, xrays);
  polarisation = cross(vectors, ott.utils.Vector(origin, ones(size(origin))));

  rays = ott.optics.geometric.Ray(vectors, power, polarisation);

  testCase.verifyEqual(size(rays), [yrays, xrays], ...
    'Size of rays incorrect');

end

function testPlot(testCase)

  xrays = 5;
  yrays = 3;

  origin = zeros(3, yrays, xrays);
  direction = rand(3, yrays, xrays);

  vectors = ott.utils.Vector(origin, direction);
  power = ones(1, yrays, xrays);
  polarisation = cross(vectors, ott.utils.Vector(origin, ones(size(origin))));
  
  rays = ott.optics.geometric.Ray(vectors, power, polarisation);
  
  h = figure();
  plot(rays);
  close(h);

end

function testScatterPlane(testCase)

  % Generate a plane to scatter rays off
  normal = [1;1;1];
  offset = sqrt(3);
  plane = ott.shapes.Plane(normal, offset);
  
  % Setup some rays
  xrays = 5;
  yrays = 3;

  origin = zeros(3, yrays, xrays);
  direction = rand(3, yrays, xrays);

  vectors = ott.utils.Vector(origin, direction);
  power = ones(1, yrays, xrays);
  polarisation = cross(vectors, ott.utils.Vector(origin, ones(size(origin))));
  
  rays = ott.optics.geometric.Ray(vectors, power, polarisation);
  
  [rrays, trays, perp] = rays.scatter(plane.normal, n1, n2);
  
  h = figure();
  plot(plane);
  hold on;
  plot(rays, 'k');
  plot(rrays, 'r');
  plot(trays, 'b');
  plot(perp, 'k');
  hold off;
  close(h);
end
  
  
