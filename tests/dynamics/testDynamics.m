function tests = testDynamics
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstructDefault(testCase)

  dynamics = ott.dynamics.Isolated();
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

  dynamics = ott.dynamics.Isolated('beam', beam, 'timeStep', timeStep, ...
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

  dynamics = ott.dynamics.Isolated('beam', beam, 'timeStep', timeStep, ...
      'temperature', temperature, 'particle', particle);

  [t, x, R] = dynamics.simulate(2*timeStep);
  testCase.verifySize(t, [1, 3], 't');
  testCase.verifySize(x, [3, 3], 'x');
  testCase.verifySize(R, [3, 9], 'R');

  dynamics.particle.mass = 1;
  [t, x, R] = dynamics.simulate(3*timeStep);
  testCase.verifySize(t, [1, 4], 'mt');
  testCase.verifySize(x, [3, 4], 'mx');
  testCase.verifySize(R, [3, 12], 'mR');

end

function testSetupAxes(testCase)

  h = figure();
  testCase.addTeardown(@() close(h));
  
  % Test no data
  sim = ott.dynamics.Isolated();
  pd = sim.setupAxes();
  testCase.verifyEqual(pd, struct('running', true, 'axes', []));
  
  ax = [];
  sim.particle.shape = ott.shape.Sphere();
  pd = sim.setupAxes(ax, 1);
  
  testCase.verifyInstanceOf(pd.axes, 'matlab.graphics.axis.Axes');
  testCase.verifyInstanceOf(pd.patch, 'matlab.graphics.primitive.Patch');
  testCase.verifyInstanceOf(pd.stopButton, 'matlab.ui.control.UIControl');
  testCase.verifyTrue(isnumeric(pd.patchVertices));
  
  % Test stop button
  testCase.verifyEqual(pd.stopButton.UserData, [], 'initial udata');
  pd.stopButton.Callback();
  testCase.verifyEqual(pd.stopButton.UserData, 'stop', 'after click');
  
  % Test existing axes
  ax = axes();
  pd = sim.setupAxes(ax, 1);
  testCase.verifyEqual(pd.axes, ax);
  
end

function testUpdatePlot(testCase)

  h = figure();
  testCase.addTeardown(@() close(h));
  
  sim = ott.dynamics.Isolated();
  pd = sim.setupAxes();
  pd2 = sim.updatePlot(pd, [], []);
  testCase.verifyEqual(pd, pd2);
  
  ax = [];
  sim.particle.shape = ott.shape.Sphere();
  pd = sim.setupAxes(ax, 1);
  pd.time = 0;
  t = 1.0;
  sim.updatePlot(pd, t, [1;0;0], eye(3));
  testCase.verifyEqual(pd.patch.Vertices, pd.patchVertices.' + [1,0,0]);

end
