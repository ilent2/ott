function tests = testVelocityField
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testOutputSizes(testCase)

  xy = randn(2, 100);
  t = linspace(0, 1, 100);

  [velocity, xedges, yedges, counts] = ott.tools.velocityField(t, xy);

  testCase.verifySize(velocity, [numel(xedges)-1, numel(yedges)-1, 2], 'velocity');
  testCase.verifySize(counts, [numel(xedges)-1, numel(yedges)-1], 'counts');

end

