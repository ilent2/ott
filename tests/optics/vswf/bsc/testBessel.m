function tests = testBscBessel
  tests = functiontests(localfunctions);
end

function setup(testCase)

  addpath('../');

end

function testConstruct(testCase)

  beam = ott.BscBessel(12, [0.0, pi/4]);

  import matlab.unittest.constraints.IsEqualTo;
  testCase.verifyThat(beam.Nbeams, IsEqualTo(2), ...
    'Incorrect number of beams stored');

end

function testMultiple(testCase)
% Combined beams should produce the same output as individual beams

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-14;
  Nmax = 10;
  
  theta1 = 0.2*pi/2;
  theta2 = 0.6*pi/2;
  theta3 = 0.5*pi/2;

  beam1 = ott.BscBessel(Nmax, theta1, 'Etheta', 1, 'Ephi', 0);
  beam2 = ott.BscBessel(Nmax, theta2, 'Etheta', 1, 'Ephi', 0);
  beam3 = ott.BscBessel(Nmax, theta3, 'Etheta', 1, 'Ephi', 0);
  beamc = ott.BscBessel(Nmax, [theta1, theta2, theta3], 'Etheta', 1, 'Ephi', 0);
  
  
  testCase.verifyThat(beamc.beam(1).getCoefficients, IsEqualTo(beam1.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'First beam doesnt match');
  
  testCase.verifyThat(beamc.beam(2).getCoefficients, IsEqualTo(beam2.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Second beam doesnt match');
  
  testCase.verifyThat(beamc.beam(3).getCoefficients, IsEqualTo(beam3.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Second beam doesnt match');

end

function testZeroAngle(testCase)
% Our toolbox will use the convention to match the output of BscPlane
% when the bessel beam goes to zero angle for linear polarisation.
%
% For spherical polarisation we need to work out what we want

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-14;
  Nmax = 10;

  % Test x polarisation
  beam = ott.BscBessel(Nmax, 0.0, 'polarisation', [1, 0]);
  beamPlane = ott.BscPlane(0, 0, 'polarisation', [1, 0], 'Nmax', Nmax);
  testCase.verifyThat(beam.getCoefficients, IsEqualTo(beamPlane.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Incorrect coefficients for x polarisation');
  
  % Test y polarisation
  beam = ott.BscBessel(Nmax, 0.0, 'polarisation', [0, 1]);
  beamPlane = ott.BscPlane(0, 0, 'polarisation', [0, 1], 'Nmax', Nmax);
  testCase.verifyThat(beam.getCoefficients, IsEqualTo(beamPlane.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Incorrect coefficients for y polarisation');
  
  % We need to handle these cases differently
  % These are non-sensical, so we won't worry about them for now
  beam = ott.BscBessel(Nmax, 0.0, 'Etheta', 1, 'Ephi', 0);
  beam = ott.BscBessel(Nmax, 0.0, 'Etheta', 0, 'Ephi', 1);

end

function testValues(testCase)
% We should also test values from the previous version of the toolbox

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-4;
  
  beam1 = ott.BscBessel(2, pi/4, 'Etheta', 1, 'Ephi', 0);
  beam2 = ott.BscBessel(2, 3*pi/4, 'Etheta', 1, 'Ephi', 0);
  
  t1a = sparse([], [], [], 8, 1);
  t2a = sparse([], [], [], 8, 1);
  t1b = sparse([2, 6], [1, 1], [3.0700i, 4.8541], 8, 1);
  t2b = sparse([2, 6], [1, 1], [3.0700i, -4.8541], 8, 1);
  
  testCase.verifyThat(beam1.getCoefficients, IsEqualTo([t1a; t1b], ...
    'Within', AbsoluteTolerance(tol)), ...
    'Incorrect value test 1');
  
  testCase.verifyThat(beam2.getCoefficients, IsEqualTo([t2a; t2b], ...
    'Within', AbsoluteTolerance(tol)), ...
    'Incorrect value test 1');
  
end