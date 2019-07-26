function tests = testFindTraps()
  % Unit tests for ott.find_traps
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  % Ensure the ott package is in our path
  addpath('../');
  
  % Disable warnings during testing
  ott.change_warnings('off');
end

function testOneSinusoid(testCase)
  % Test a ideal sinusoid trap

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  x = linspace(-1, 1);
  x(1) = []; x(end) = [];  % Drop zeros at ends
  y = -2*sin(pi*x);

  traps = ott.find_traps(x, y);
  
  tol1 = 1e-2;
  dy = @(x) -2*pi*cos(pi*x);

  testCase.verifyThat(numel(traps), IsEqualTo(1), ...
      'Wrong number of traps identified');
  testCase.verifyThat(traps(1).stiffness, IsEqualTo(dy(0), ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap stiffness incorrect');
  testCase.verifyThat(traps(1).position, IsEqualTo(0, ...
      'Within', AbsoluteTolerance(1e-6)), ...
      'Trap position incorrect');
  testCase.verifyThat(traps(1).depth, IsEqualTo(2, ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap depth incorrect');
  testCase.verifyThat(traps(1).minmax_force, IsEqualTo([2, -2], ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap minmax_force incorrect');
  testCase.verifyThat(traps(1).minmax_position, IsEqualTo([-0.5, 0.5], ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap minmax_position incorrect');
  testCase.verifyThat(traps(1).range, IsEqualTo([-Inf, Inf], ...
      'Within', AbsoluteTolerance(1e-6)), ...
      'Trap range incorrect');
end

function testTwoSinusoid(testCase)
  % Test a ideal sinusoid trap

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  x = linspace(-2, 2, 200);
  x(1) = []; x(end) = [];  % Drop zeros at ends
  y = 2*sin(pi*x);

  traps = ott.find_traps(x, y);
  
  tol1 = 1e-2;
  dy = @(x) 2*pi*cos(pi*x);

  testCase.verifyThat(numel(traps), IsEqualTo(2), ...
      'Wrong number of traps identified');

  testCase.verifyThat(traps(1).stiffness, IsEqualTo(dy(-1), ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap1 stiffness incorrect');
  testCase.verifyThat(traps(2).stiffness, IsEqualTo(dy(1), ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap2 stiffness incorrect');

  testCase.verifyThat(traps(1).position, IsEqualTo(-1, ...
      'Within', AbsoluteTolerance(1e-6)), ...
      'Trap1 position incorrect');
  testCase.verifyThat(traps(2).position, IsEqualTo(1, ...
      'Within', AbsoluteTolerance(1e-6)), ...
      'Trap2 position incorrect');

  testCase.verifyThat(traps(1).depth, IsEqualTo(2, ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap1 depth incorrect');
  testCase.verifyThat(traps(2).depth, IsEqualTo(2, ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap2 depth incorrect');

  testCase.verifyThat(traps(1).minmax_force, IsEqualTo([2, -2], ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap1 minmax_force incorrect');
  testCase.verifyThat(traps(2).minmax_force, IsEqualTo([2, -2], ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap2 minmax_force incorrect');

  testCase.verifyThat(traps(1).minmax_position, IsEqualTo([-1.5, -0.5], ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap1 minmax_position incorrect');
  testCase.verifyThat(traps(2).minmax_position, IsEqualTo([0.5, 1.5], ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap2 minmax_position incorrect');

  testCase.verifyThat(traps(1).range, IsEqualTo([-Inf, 0], ...
      'Within', AbsoluteTolerance(1e-6)), ...
      'Trap1 range incorrect');
  testCase.verifyThat(traps(2).range, IsEqualTo([0, Inf], ...
      'Within', AbsoluteTolerance(1e-6)), ...
      'Trap2 range incorrect');
end

function testOneNegativeSinusoid(testCase)
  % Test a ideal sinusoid trap

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  x = linspace(-1, 1);
  x(1) = []; x(end) = [];  % Drop zeros at ends
  y = 2*sin(pi*x);

  traps = ott.find_traps(x, y, 'keep_unstable', true);

  tol1 = 1e-2;
  dy = @(x) 2*pi*cos(pi*x);

  testCase.verifyThat(numel(traps), IsEqualTo(1), ...
      'Wrong number of traps identified');
  testCase.verifyThat(traps(1).stiffness, IsEqualTo(dy(0), ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap stiffness incorrect');
  testCase.verifyThat(traps(1).position, IsEqualTo(0, ...
      'Within', AbsoluteTolerance(1e-6)), ...
      'Trap position incorrect');
  testCase.verifyThat(traps(1).depth, IsEqualTo(2, ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap depth incorrect');
  testCase.verifyThat(traps(1).minmax_force, IsEqualTo([-2, 2], ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap minmax_force incorrect');
  testCase.verifyThat(traps(1).minmax_position, IsEqualTo([-0.5, 0.5], ...
      'Within', AbsoluteTolerance(tol1)), ...
      'Trap minmax_position incorrect');
  testCase.verifyThat(traps(1).range, IsEqualTo([-Inf, Inf], ...
      'Within', AbsoluteTolerance(1e-6)), ...
      'Trap range incorrect');
end

function testSawTooth(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1e-2;

  x = [0, 1, 2, 3, 4];
  y = [1, -1, 1, -1, 1];
  
  % Test all equilibrium 
  
  eqs = [0.5, 1.5, 2.5, 3.5];
  
  traps = ott.find_traps(x, y, 'keep_unstable', true);
  
  testCase.verifyThat(numel(traps), IsEqualTo(length(eqs)), ...
      'Wrong number of stable+unstable traps identified');
  testCase.verifyThat([traps(2:3).position], IsEqualTo(eqs(2:3), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Trap position incorrect (keep_unstable)');
  
  % Test stable equilibria
  
  eqs = [0.5, 2.5];

  traps = ott.find_traps(x, y, 'keep_unstable', false);
  
  testCase.verifyThat(numel(traps), IsEqualTo(length(eqs)), ...
      'Wrong number of stable traps identified');
  testCase.verifyThat([traps(2).position], IsEqualTo(eqs(2), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Trap position incorrect');
end

function testFlatCentreTrap(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1e-2;
  
  y = [0.1, 0.5, 0, 0, -0.5, -0.1];
  x = [-3, -2, -1, 1, 2, 3];
  
  eqs = [0.0];
  
  traps = ott.find_traps(x, y, 'keep_unstable', true);
  
  testCase.verifyThat(numel(traps), IsEqualTo(length(eqs)), ...
      'Wrong number of traps identified with keep_unstable=true');
  
  % Test stable equilibria

  traps = ott.find_traps(x, y, 'keep_unstable', false);
  
  testCase.verifyThat(numel(traps), IsEqualTo(length(eqs)), ...
      'Wrong number of traps identified with keep_unstable=false');
  testCase.verifyThat([traps.position], IsEqualTo(eqs, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Trap position incorrect');
    
% Do we want to have two definitions of stiffness?
%   testCase.verifyThat([traps.stiffness], IsEqualTo(0, ...
%       'Within', AbsoluteTolerance(tol)), ...
%       'Trap stiffness incorrect');

end

function testFlatCentrePoly3Trap(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1e-2;
  
  x = linspace(-1, 1, 100);
  p = [-1, 0, 0, 0];
  y = polyval(p, x);
  
  eqs = [0.0];
  
  traps = ott.find_traps(x, y, 'keep_unstable', true);
  
  testCase.verifyThat(numel(traps), IsEqualTo(length(eqs)), ...
      'Wrong number of traps identified with keep_unstable=true');
  
  % Test stable equilibria

  traps = ott.find_traps(x, y, 'keep_unstable', false);
  
  testCase.verifyThat(numel(traps), IsEqualTo(length(eqs)), ...
      'Wrong number of traps identified with keep_unstable=false');
  testCase.verifyThat([traps.position], IsEqualTo(eqs, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Trap position incorrect');
  testCase.verifyThat([traps.stiffness], IsEqualTo(0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Trap position incorrect');

end

function testSingleTrapGroup(testCase)

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
