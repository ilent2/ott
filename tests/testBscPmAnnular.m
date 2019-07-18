function tests = testBscPmAnnular
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../');
end

function testUniform(testCase)
  
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 0.1;

  beam = ott.BscPmAnnular([0.5, 1.0]);
  
  beam.basis = 'incoming';
  E = beam.farfield([pi, pi-0.8*pi/2], [0, 0]);
  E = sqrt(sum(abs(E).^2, 1));
  
  testCase.verifyThat(E(1), ...
    IsEqualTo(0.0, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Centre of beam does not go to zero');
  
  testCase.verifyThat(E(2), ...
    IsEqualTo(1.0, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Edge of beam does not go to one');
  
  figure();
  beam.basis = 'incoming';
  beam.visualiseFarfield('dir', 'neg');

end

function testFromProfile(testCase)

  profile = linspace(0, 1, 20);
  beam = ott.BscPmAnnular([0.5, 1.0], 'profile', profile);
  
  figure();
  beam.basis = 'incoming';
  beam.visualiseFarfield('dir', 'neg');
  
end

function testFromBeam(testCase)

  oldBeam = ott.BscPmGauss();
  oldBeam.basis = 'incoming';

  beam = ott.BscPmAnnular([0.5, 1.0], 'profile', oldBeam);
  
  figure();
  beam.basis = 'incoming';
  beam.visualiseFarfield('dir', 'neg');
  
end
