function tests = testBesselBasis
  tests = functiontests(localfunctions);
end

function setup(testCase)
  addpath('../../../');
end

% See testConsistency for default constructor

function testConstructOptional(testCase)

  Nmax = 1;
  angle = 0.1;

  beam = ott.beam.vswf.BesselBasis(Nmax, angle);

  testCase.verifyEqual(beam.Nmax, Nmax, 'nmax');
  testCase.verifyEqual(beam.angle, angle, 'angle');

end

% TODO: review tests from here down

function testZeroAngle(testCase)
% Our toolbox will use the convention to match the output of BscPlane
% when the bessel beam goes to zero angle for linear polarisation.
%
% For spherical polarisation we need to work out what we want
%
% TODO: Do we want to change this in version 2?

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  import ott.optics.vswf.*;
  tol = 1.0e-14;
  Nmax = 10;

  % Test x polarisation
  beam = bsc.Bessel(Nmax, 0.0, 'polarisation', [1, 0]);
  beamPlane = bsc.Plane(0, 0, 'polarisation', [1, 0], 'Nmax', Nmax);
  testCase.verifyThat(beam.getCoefficients, IsEqualTo(beamPlane.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Incorrect coefficients for x polarisation');
  
  % Test y polarisation
  beam = bsc.Bessel(Nmax, 0.0, 'polarisation', [0, 1]);
  beamPlane = bsc.Plane(0, 0, 'polarisation', [0, 1], 'Nmax', Nmax);
  testCase.verifyThat(beam.getCoefficients, IsEqualTo(beamPlane.getCoefficients, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Incorrect coefficients for y polarisation');
  
  % We need to handle these cases differently
  % These are non-sensical, so we won't worry about them for now
  beam = bsc.Bessel(Nmax, 0.0, 'Etheta', 1, 'Ephi', 0);
  beam = bsc.Bessel(Nmax, 0.0, 'Etheta', 0, 'Ephi', 1);

end

function testValues(testCase)
% We should also test values from the previous version of the toolbox

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  import ott.optics.vswf.*;
  tol = 1.0e-4;
  
  beam1 = bsc.Bessel(2, pi/4, 'Etheta', 1, 'Ephi', 0);
  beam2 = bsc.Bessel(2, 3*pi/4, 'Etheta', 1, 'Ephi', 0);
  
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
