function tests = testDynamics
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructDefault(testCase)

  dynamics = ott.tools.Dynamics();
  testCase.verifyEqual(dynamics.particle, ott.particle.Fixed(), 'particle');
  testCase.verifyEqual(dynamics.timeStep, 1e-4, 'time');
  testCase.verifyEqual(dynamics.beam, ott.beam.Empty(), 'beam');
  testCase.verifyEqual(dynamics.temperature, 300, 'temperature');

end

function testConstructArgs(testCase)

  beam = ott.beam.Gaussian();
  particle = ott.particle.Variable();
  timeStep = 1;
  temperature = 1;

  dynamics = ott.tools.Dynamics('beam', beam, 'timeStep', timeStep, ...
      'temperature', temperature, 'particle', particle);
  testCase.verifyEqual(dynamics.particle, particle, 'particle');
  testCase.verifyEqual(dynamics.timeStep, timeStep, 'time');
  testCase.verifyEqual(dynamics.beam, beam, 'beam');
  testCase.verifyEqual(dynamics.temperature, temperature, 'temperature');

end

function testSimulateCoverage(testCase)

  beam = ott.beam.Gaussian();
  particle = ott.particle.Fixed();
  timeStep = 1;
  temperature = 1;

  dynamics = ott.tools.Dynamics('beam', beam, 'timeStep', timeStep, ...
      'temperature', temperature, 'particle', particle);

  [t, x, R] = dynamics.simulate(2*timeStep);
  testCase.verifySize(t, [1, 2], 't');
  testCase.verifySize(x, [3, 2], 'x');
  testCase.verifySize(R, [3, 6], 'R');

  dynamics.particle.mass = 1;
  [t, x, R] = dynamics.simulate(3*timeStep);
  testCase.verifySize(t, [1, 3], 'mt');
  testCase.verifySize(x, [3, 3], 'mx');
  testCase.verifySize(R, [3, 9], 'mR');

end

