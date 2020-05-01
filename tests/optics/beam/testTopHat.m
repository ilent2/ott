function tests = testTopHat
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstruct(testCase)

  profile = ott.optics.beam.TopHat.ProfileSquare(1.0, [1;0]);
  direction = randn(3, 1);
  polarisation = randn(3, 1);
  field = 1.0;
  origin = [0;0;0];
  power = 1.0;
  beam = ott.optics.beam.TopHat('profile', profile, ...
    'direction', direction, 'field', field, 'origin', origin, ...
    'polarisation', polarisation, 'power', power);
  
  testCase.verifyEqual(beam.profile, profile, 'profile');
  testCase.verifyEqual(beam.direction, direction, 'direction');
  testCase.verifyEqual(beam.polarisation, polarisation, 'polarisation');
  testCase.verifyEqual(beam.field, field, 'field');
  testCase.verifyEqual(beam.origin, origin, 'origin');
  testCase.verifyEqual(beam.power, power, 'power');

end

function testPowerEstimationSymbolic(testCase)

  field = 1.5;
  width = 2.0;
  profile = ott.optics.beam.TopHat.ProfileSquare(width, [1;0]);
  targetPower = field * width.^2;
  
  beam = ott.optics.beam.TopHat('profile', profile, ...
    'field', field, 'power', 'symbolic');
  
  testCase.verifyEqual(beam.power, targetPower, 'power');

end

function testPowerEstimationNumeric(testCase)

  field = 1.5;
  width = 2.0;
  profile = ott.optics.beam.TopHat.ProfileSquare(width, [1;0]);
  targetPower = field * width.^2;
  
  beam = ott.optics.beam.TopHat('profile', profile, ...
    'field', field, 'power', 'numeric');
  
  testCase.verifyEqual(beam.power, targetPower, 'AbsTol', 1.0e-2, 'power');

end

function a = evaluateProfileArea(profile)

% Doesn't work nicely for circles...
%   func = piecewise(sym(profile), 1, 0);
%   a = double(int(int(func, -Inf, Inf), -Inf, Inf));

  a = integral2(@(x, y) double(profile(x, y)), -Inf, Inf, -Inf, Inf);

end

function testProfileSquare(testCase)

  width = 1.0;
  fcn = ott.optics.beam.TopHat.ProfileSquare(width, [1; 0]);
  
  xrange = linspace(-2, 2, 50);
  yrange = linspace(-2, 2, 50);
  [xx, yy] = meshgrid(xrange, yrange);
  
  b = fcn(xx(:).', yy(:).');
  target = abs(xx) <= width/2 & abs(yy) <= width/2;
  
%   figure();
%   imagesc(reshape(b, size(xx)));
%   axis image;
  
  testCase.verifyEqual(b, target(:).', 'val');
  testCase.verifyEqual(evaluateProfileArea(fcn), width.^2, ...
    'AbsTol', 1.0e-2, 'area');
  
  width = 1.0;
  fcn = ott.optics.beam.TopHat.ProfileSquare(width, [1; 1]);
  
  xy = [1.1, 1.1, -1.1, -1.1, 0.0; 1.1, -1.1, 1.1, -1.1, 0.0];
  target = [false, false, false, false, true];
  
  b = fcn(xy(1, :), xy(2, :));
  testCase.verifyEqual(b, target, 'diag');
  testCase.verifyEqual(evaluateProfileArea(fcn), width.^2, ...
    'AbsTol', 1.0e-2, 'diag area');
  
%   figure();
%   imagesc(reshape(fcn(xx(:).', yy(:).'), size(xx)));
%   axis image;

end

function testProfileCircle(testCase)

  radius = 1.0;
  fcn = ott.optics.beam.TopHat.ProfileCircle(radius);
  
  xrange = linspace(-2, 2, 50);
  yrange = linspace(-2, 2, 50);
  [xx, yy] = meshgrid(xrange, yrange);
  
  b = fcn(xx(:).', yy(:).');
  target = xx.^2 + yy.^2 <= radius.^2;
  
  testCase.verifyEqual(b, target(:).', 'val');
  testCase.verifyEqual(evaluateProfileArea(fcn), pi*radius.^2, ...
    'AbsTol', 1.0e-2, 'area');
  
%   figure();
%   imagesc(reshape(b, size(xx)));
%   axis image;

end

function testProfileEllipse(testCase)

  radii = [0.5, 1.0];
  fcn = ott.optics.beam.TopHat.ProfileEllipse(radii, [1; 0]);
  
  xrange = linspace(-2, 2, 50);
  yrange = linspace(-2, 2, 50);
  [xx, yy] = meshgrid(xrange, yrange);
  
  b = fcn(xx(:).', yy(:).');
  
  testCase.verifyEqual(size(b), [1, numel(xx)], 'val');
  testCase.verifyEqual(evaluateProfileArea(fcn), pi*prod(radii), ...
    'AbsTol', 1.0e-2, 'area');
  
%   figure();
%   imagesc(reshape(b, size(xx)));
%   axis image;

end

function testProfileRectangle(testCase)

  widths = [0.5, 1.0];
  fcn = ott.optics.beam.TopHat.ProfileRectangle(widths, [1; 0]);
  
  xrange = linspace(-2, 2, 50);
  yrange = linspace(-2, 2, 50);
  [xx, yy] = meshgrid(xrange, yrange);
  
  b = fcn(xx(:).', yy(:).');
  
  testCase.verifyEqual(size(b), [1, numel(xx)], 'val');
  testCase.verifyEqual(evaluateProfileArea(fcn), prod(widths), ...
    'AbsTol', 1.0e-2, 'area');
  
%   figure();
%   imagesc(reshape(b, size(xx)));
%   axis image;

end

function testCastPlane(testCase)

  beam = ott.optics.beam.TopHat();
  beam2 = ott.optics.beam.PlaneWave(beam);

end

function testVisualise(testCase)

  field = @(x, y) exp(-x.^2 - y.^2);
  beam = ott.optics.beam.TopHat('field', field);
  
  h = figure();
  beam.visualise();
  close(h);

end

