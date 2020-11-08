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
  y = -2*sin(3*pi*x);

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', true);

  h1 = plot(x, y);
  hold on;
  h2 = traps.plot();
  hold off;

  testCase.verifyTrue(all(ishandle(h1)), 'h1');
  testCase.verifyTrue(all(ishandle(h2)), 'h2');

end

function testBisection(testCase)

  x = [-0.5, 0.5];
  y = @(x) -2*sin(pi*x);

  traps = ott.tools.FindTraps1d.Bisection(x, y);
  
  testCase.verifyEqual(traps.position, 0, 'p');
  testCase.verifyEqual(traps.stiffness, -2*pi, 's');

end

function testNewton(testCase)

  x = 0.2;
  y = @(x) -2*sin(pi*x);

  traps = ott.tools.FindTraps1d.Newton(x, y);
  
  testCase.verifyWarningFree(@() ott.tools.FindTraps1d.Newton(x, y));
  
  testCase.verifyEqual(traps.position, 0, 'AbsTol', 1e-15, 'p');
  testCase.verifyEqual(traps.stiffness, -2*pi, 's');

end

function testSingleTrapGroup(testCase)

  y = [0.1, 1, -1, -0.1];
  x = 1:numel(y);

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', true);

  testCase.verifySize(traps, [1, 1], ...
      'Wrong number of traps identified with keep_unstable=true');

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', false);

  testCase.verifySize(traps, [1, 1], ...
      'Wrong number of traps identified with keep_unstable=false');

  [~, stats] = traps.groupStable();

  testCase.verifySize(numel(stats), [1,1], ...
      'Wrong number of traps identified with group_stable=true');

end

function testDoubleTrapGroup(testCase)

  y = [0.1, 1, -1, 0.5, -0.5, 0.3, -2, -1];
  x = 1:numel(y);

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', true);

  testCase.verifySize(traps, [1, 5], ...
      'Wrong number of traps identified with keep_unstable=true');

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', false);

  testCase.verifySize(traps, [1, 3], ...
      'Wrong number of traps identified with keep_unstable=false');

  [~, stats] = traps.groupStable();

  testCase.verifySize(numel(stats), [1,1], ...
      'Wrong number of traps identified with group_stable=true');

end

