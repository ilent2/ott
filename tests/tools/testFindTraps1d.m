function tests = testFindTraps1d
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  traps = ott.tools.FindTraps1d;

end

function testPlot(testCase)

  f = figure();
  testCase.addTeardown(@close, f);

  x = linspace(-1, 1);
  x(1) = []; x(end) = [];  % Drop zeros at ends
  y = -2*sin(pi*x);

  traps = ott.tools.FindTraps1d.FromArray(x, y);

  h1 = plot(x, y);
  hold on;
  h2 = traps.plot();
  hold off;

  testCase.verifyThat(ishandle(h1), 'h1');
  testCase.verifyThat(ishandle(h2), 'h2');

end

function testSingleTrapGroup(testCase)

  % TODO:Rewrite

  import matlab.unittest.constraints.IsEqualTo;

  y = [0.1, 1, -1, -0.1];
  x = 1:numel(y);

  traps = ott.find_traps(x, y, 'keep_unstable', true);

  testCase.verifyThat(numel(traps), IsEqualTo(1), ...
      'Wrong number of traps identified with keep_unstable=true');

  traps = ott.find_traps(x, y, 'keep_unstable', false);

  testCase.verifyThat(numel(traps), IsEqualTo(1), ...
      'Wrong number of traps identified with keep_unstable=false');

  traps = ott.find_traps(x, y, 'keep_unstable', false, 'group_stable', true);

  testCase.verifyThat(numel(traps), IsEqualTo(1), ...
      'Wrong number of traps identified with group_stable=true');

end

function testDoubleTrapGroup(testCase)

  % TODO:Rewrite

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1e-2;

  y = [0.1, 1, -1, 0.5, -0.5, 0.3, -2, -1];
  x = 1:numel(y);

  traps = ott.find_traps(x, y, 'keep_unstable', true);

  testCase.verifyThat(numel(traps), IsEqualTo(5), ...
      'Wrong number of traps identified with keep_unstable=true');

  traps = ott.find_traps(x, y, 'keep_unstable', false);

  testCase.verifyThat(numel(traps), IsEqualTo(3), ...
      'Wrong number of traps identified with keep_unstable=false');

  traps = ott.find_traps(x, y, 'keep_unstable', false, 'group_stable', true);

  testCase.verifyThat(numel(traps), IsEqualTo(1), ...
      'Wrong number of traps identified with group_stable=true');

end

