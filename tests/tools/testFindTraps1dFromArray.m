function tests = testFindTraps1dFromArray
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testOneSinusoid(testCase)
  % Test a ideal sinusoid trap

  x = linspace(-1, 1);
  x(1) = []; x(end) = [];  % Drop zeros at ends
  y = -2*sin(pi*x);

  traps = ott.tools.FindTraps1d.FromArray(x, y);

  dy = @(x) -2*pi*cos(pi*x);

  testCase.assertSize(traps, [1, 1], ...
      'Wrong number of traps identified');
  testCase.verifyEqual(traps(1).stiffness, dy(0), ...
      'AbsTol', 1e-2, 'Trap stiffness incorrect');
  testCase.verifyEqual(traps(1).position, 0, ...
      'AbsTol', 1e-6, 'Trap position incorrect');
  testCase.verifyEqual(traps(1).depth, 2, ...
      'AbsTol', 1e-2, 'Trap depth incorrect');
  testCase.verifyEqual(traps(1).minforce, 2, ...
      'AbsTol', 1e-2, 'Trap minmax_force incorrect');
  testCase.verifyEqual(traps(1).maxforce, -2, ...
      'AbsTol', 1e-2, 'Trap minmax_force incorrect');
  testCase.verifyEqual(traps(1).minposition, -0.5, ...
      'AbsTol', 1e-2, 'Trap minmax_position incorrect');
  testCase.verifyEqual(traps(1).maxposition, 0.5, ...
      'AbsTol', 1e-2, 'Trap minmax_position incorrect');
  testCase.verifyEqual(traps(1).range, [-Inf, Inf], ...
      'AbsTol', 1e-6, 'Trap range incorrect');
end

function testTwoSinusoid(testCase)
  % Test a ideal sinusoid trap

  x = linspace(-2, 2, 200);
  x(1) = []; x(end) = [];  % Drop zeros at ends
  y = 2*sin(pi*x);

  traps = ott.tools.FindTraps1d.FromArray(x, y);
  dy = @(x) 2*pi*cos(pi*x);

  testCase.assertSize(traps, [1, 2], ...
      'Wrong number of traps identified');

  testCase.verifyEqual(traps(1).stiffness, dy(-1), ...
      'AbsTol', 1e-2, 'Trap1 stiffness incorrect');
  testCase.verifyEqual(traps(2).stiffness, dy(1), ...
      'AbsTol', 1e-2, 'Trap2 stiffness incorrect');

  testCase.verifyEqual(traps(1).position, -1, ...
      'AbsTol', 1e-6, 'Trap1 position incorrect');
  testCase.verifyEqual(traps(2).position, 1, ...
      'AbsTol', 1e-6, 'Trap2 position incorrect');

  testCase.verifyEqual(traps(1).depth, 2, ...
      'AbsTol', 1e-2, 'Trap1 depth incorrect');
  testCase.verifyEqual(traps(2).depth, 2, ...
      'AbsTol', 1e-2, 'Trap2 depth incorrect');

  testCase.verifyEqual([traps.minforce], [2, 2], ...
      'AbsTol', 1e-2, 'Trap1 minmax_force incorrect');
  testCase.verifyEqual([traps.maxforce], [-2, -2], ...
      'AbsTol', 1e-2, 'Trap2 minmax_force incorrect');

  testCase.verifyEqual([traps.minposition], [-1.5, 0.5], ...
      'AbsTol', 1e-2, 'Trap1 minmax_position incorrect');
  testCase.verifyEqual([traps.maxposition], [-0.5, 1.5], ...
      'AbsTol', 1e-2, 'Trap2 minmax_position incorrect');

  testCase.verifyEqual(traps(1).range, [-Inf, 0], ...
      'AbsTol', 1e-6, 'Trap1 range incorrect');
  testCase.verifyEqual(traps(2).range, [0, Inf], ...
      'AbsTol', 1e-6, 'Trap2 range incorrect');
end

function testOneNegativeSinusoid(testCase)
  % Test a ideal sinusoid trap

  x = linspace(-1, 1);
  x(1) = []; x(end) = [];  % Drop zeros at ends
  y = 2*sin(pi*x);

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', true);

  dy = @(x) 2*pi*cos(pi*x);

  testCase.assertSize(traps, [1, 1], ...
      'Wrong number of traps identified');
  testCase.verifyEqual(traps(1).stiffness, dy(0), ...
      'AbsTol', 1e-2, 'Trap stiffness incorrect');
  testCase.verifyEqual(traps(1).position, 0, ...
      'AbsTol', 1e-2, 'Trap position incorrect');
  testCase.verifyEqual(traps(1).depth, 2, ...
      'AbsTol', 1e-2, 'Trap depth incorrect');
  testCase.verifyEqual(traps(1).minforce, -2, ...
      'AbsTol', 1e-2, 'Trap minmax_force incorrect');
  testCase.verifyEqual(traps(1).maxforce, 2, ...
      'AbsTol', 1e-2, 'Trap minmax_force incorrect');
  testCase.verifyEqual(traps(1).minposition, -0.5, ...
      'AbsTol', 1e-2, 'Trap minmax_position incorrect');
  testCase.verifyEqual(traps(1).maxposition, 0.5, ...
      'AbsTol', 1e-2, 'Trap minmax_position incorrect');
  testCase.verifyEqual(traps(1).range, [-Inf, Inf], ...
      'AbsTol', 1e-2, 'Trap range incorrect');
end

function testSawTooth(testCase)

  x = [0, 1, 2, 3, 4];
  y = [1, -1, 1, -1, 1];

  % Test all equilibrium

  eqs = [0.5, 1.5, 2.5, 3.5];

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', true);

  testCase.assertSize(traps, size(eqs), ...
      'Wrong number of stable+unstable traps identified');
  testCase.verifyEqual([traps(2:3).position], eqs(2:3), ...
      'AbsTol', 1e-2, 'Trap position incorrect (keep_unstable)');

  % Test stable equilibria

  eqs = [0.5, 2.5];

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', false);

  testCase.assertSize(traps, size(eqs), ...
      'Wrong number of stable traps identified');
  testCase.verifyEqual([traps(2).position], eqs(2), ...
      'AbsTol', 1e-2, 'Trap position incorrect');
end

function testFlatCentreTrap(testCase)

  y = [0.1, 0.5, 0, 0, -0.5, -0.1];
  x = [-3, -2, -1, 1, 2, 3];

  eqs = 0.0;

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', true);

  testCase.verifySize(traps, size(eqs), ...
      'Wrong number of traps identified with keep_unstable=true');

  % Test stable equilibria

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', false);

  testCase.verifySize(traps, size(eqs), ...
      'Wrong number of traps identified with keep_unstable=false');
  testCase.verifyEqual([traps.position], eqs, ...
      'AbsTol', 1e-2, 'Trap position incorrect');

% Do we want to have two definitions of stiffness?
%   testCase.verifyThat([traps.stiffness], IsEqualTo(0, ...
%       'Within', AbsoluteTolerance(tol)), ...
%       'Trap stiffness incorrect');

end

function testFlatCentrePoly3Trap(testCase)

  x = linspace(-1, 1, 100);
  p = [-1, 0, 0, 0];
  y = polyval(p, x);

  eqs = 0.0;

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', true);

  testCase.verifySize(traps, size(eqs), ...
      'Wrong number of traps identified with keep_unstable=true');

  % Test stable equilibria

  traps = ott.tools.FindTraps1d.FromArray(x, y, 'keep_unstable', false);

  testCase.assertSize(traps, size(eqs), ...
      'Wrong number of traps identified with keep_unstable=false');
  testCase.verifyEqual([traps.position], eqs, ...
      'AbsTol', 1e-2, 'Trap position incorrect');
  testCase.verifyEqual([traps.stiffness], 0, ...
      'Abstol', 1e-2, 'Trap stiffness incorrect');

end

