function tests = testBscPmAnnular
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testUniform(testCase)
  
  import ott.optics.vswf.*;
  
  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 0.1;

  beam = bsc.PmAnnular([0.5, 1.0]);
  
  beam.basis = 'incoming';
  E = beam.efarfield([[1, 1].', [pi, pi-0.8*pi/2].', [0, 0].']);
  E = sqrt(sum(abs(E).^2, 1));
  
  testCase.verifyThat(E(1), ...
    IsEqualTo(0.0, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Centre of beam does not go to zero');
  
  testCase.verifyThat(E(2), ...
    IsEqualTo(1.0, ...
    'Within', AbsoluteTolerance(tol)), ...
    'Edge of beam does not go to one');
  
  h = figure();
  beam.basis = 'incoming';
  beam.visualiseFarfield('dir', 'neg');
  close(h);

end

function testFromProfile(testCase)

  import ott.optics.vswf.*;
  
  profile = linspace(0, 1, 20);
  beam = bsc.PmAnnular([0.5, 1.0], 'profile', profile);
  
  h = figure();
  beam.basis = 'incoming';
  beam.visualiseFarfield('dir', 'neg');
  close(h);
  
end

function testFromBeam(testCase)

  import ott.optics.vswf.*;

  oldBeam = bsc.PmGauss();
  oldBeam.basis = 'incoming';

  beam = bsc.PmAnnular([0.5, 1.0], 'profile', oldBeam);
  
  h = figure();
  beam.basis = 'incoming';
  beam.visualiseFarfield('dir', 'neg');
  close(h);
  
end
