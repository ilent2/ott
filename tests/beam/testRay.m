function tests = testRay
  % Test cases based on example_ray.m script from OTGO
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function test1dArray(testCase)

  num_rays = 5;
  origin = zeros(3, num_rays);
  direction = rand(3, num_rays);

  vectors = ott.utils.Vector(origin, direction);
  power = ones(1, num_rays);
  polarisation = cross(vectors, ott.utils.Vector(origin, ones(size(origin))));

  rays = ott.beam.Ray('vector', vectors, 'field', power, 'polarisation', polarisation);

  testCase.verifyEqual(size(rays), [1, num_rays], ...
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

  rays = ott.beam.Ray('vector', vectors, 'field', power, 'polarisation', polarisation);

  testCase.verifyEqual(size(rays), [1, yrays, xrays], ...
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
  
  rays = ott.beam.Ray('vector', vectors, 'field', power, 'polarisation', polarisation);
  
  h = figure();
  plot(rays);
  close(h);

end

function testScatterPlane(testCase)

  % Generate a plane to scatter rays off
  normal = -[1;1;1];
  offset = -sqrt(3);
  index_relative = 1.33;
  plane = ott.scat.geometric.Plane(normal, index_relative, ...
      'offset', offset);
  
  % Setup some rays
  xrays = 5;
  yrays = 3;

  origin = zeros(3, yrays, xrays);
  direction = rand(3, yrays, xrays);

  vectors = ott.utils.Vector(origin, direction);
  power = ones(1, yrays, xrays);
  polarisation = cross(vectors, ott.utils.Vector(origin, ones(size(origin))));
  
  rays = ott.beam.Ray('vector', vectors, ...
    'field', power, 'polarisation', polarisation);
  
  [rrays, trays] = plane.scatter(rays);
  
  h = figure();
  plane.surf();
  hold on;
  plot(rays, 'Color', 'k');
  plot(rrays, 'Color', 'r');
  plot(trays, 'Color', 'b');
  hold off;
  close(h);
end
  
  
