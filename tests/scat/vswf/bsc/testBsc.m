function tests = bsc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testConstructDefault(testCase)

  beam = ott.optics.vswf.bsc.Bsc();
  
  testCase.verifyEqual(beam.a, [], 'a');
  testCase.verifyEqual(beam.b, [], 'b');
  testCase.verifyEqual(beam.basis, 'regular', 'basis');
  testCase.verifyEqual(beam.wavenumber, 2*pi, 'a');
  testCase.verifyEqual(beam.omega, 2*pi, 'b');
  testCase.verifyEqual(beam.dz, 0, 'dz');

  dz = 1.0;
  omega = 2.0;
  k_medium = 3.0;
  beam = ott.optics.vswf.bsc.Bsc(...
    'dz', dz, 'omega', omega, 'k_medium', k_medium);
  
  testCase.verifyEqual(beam.a, [], 'a');
  testCase.verifyEqual(beam.b, [], 'b');
  testCase.verifyEqual(beam.basis, 'regular', 'basis');
  testCase.verifyEqual(beam.wavenumber, k_medium, 'a');
  testCase.verifyEqual(beam.omega, omega, 'b');
  testCase.verifyEqual(beam.dz, dz, 'dz');
end

function testConstructFromBsc(testCase)

  bsc = ott.optics.vswf.bsc.Bsc('dz', 1.0, 'omega', 4, 'k_medium', 1);
  
  beam = ott.optics.vswf.bsc.Bsc(bsc);
  
  testCase.verifyEqual(beam.a, bsc.a, 'a');
  testCase.verifyEqual(beam.b, bsc.b, 'b');
  testCase.verifyEqual(beam.basis, bsc.basis, 'basis');
  testCase.verifyEqual(beam.wavenumber, bsc.wavenumber, 'a');
  testCase.verifyEqual(beam.omega, bsc.omega, 'b');
  testCase.verifyEqual(beam.dz, bsc.dz, 'dz');
end

function testConstructAbBasis(testCase)

  a = [];
  b = [];
  basis = 'incoming';

  beam = ott.optics.vswf.bsc.Bsc(a, b, basis);
  
  testCase.verifyEqual(beam.a, a, 'a');
  testCase.verifyEqual(beam.b, b, 'b');
  testCase.verifyEqual(beam.basis, basis, 'basis');
  testCase.verifyEqual(beam.wavenumber, 2*pi, 'a');
  testCase.verifyEqual(beam.omega, 2*pi, 'b');
  testCase.verifyEqual(beam.dz, 0, 'dz');

  beam = ott.optics.vswf.bsc.Bsc(a, b);
  
  testCase.verifyEqual(beam.a, a, 'a');
  testCase.verifyEqual(beam.b, b, 'b');
  testCase.verifyEqual(beam.basis, 'regular', 'basis');
  testCase.verifyEqual(beam.wavenumber, 2*pi, 'a');
  testCase.verifyEqual(beam.omega, 2*pi, 'b');
  testCase.verifyEqual(beam.dz, 0, 'dz');
end

function testTranslation(testCase)
  % Check that the translation functions run

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  beam = ott.optics.vswf.bsc.PmGauss();
  beam.power = 1.0;
  dz = pi/2;
  tol = 1.0e-6;

  tbeam1 = beam.translateZ(dz);

  [~, AB] = beam.translateZ(dz);
  tbeam2 = AB * beam;

  [~, A, B] = beam.translateZ(dz);
  tbeam3 = [ A B ; B A ] * beam;
  tbeam4 = beam.translate(A, B);

  % Check all beams are equal
  target = tbeam1.getCoefficients;
  testCase.verifyThat(tbeam2.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'AB * beam does not match translateZ');
  testCase.verifyThat(tbeam3.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    '[A, B; B, A] * beam does not match translateZ');
  testCase.verifyThat(tbeam4.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'beam.translate(A, B) does not match translateZ');

  xbeam1 = beam.translateXyz([dz; dz; 0]);

  [~, AB] = beam.translateXyz([dz; dz; 0]);
  xbeam2 = AB * beam;

  [~, A, B] = beam.translateXyz([dz; dz; 0]);
  xbeam3 = beam.translate(A, B);

  [~, A, B, D] = beam.translateXyz([dz; dz; 0]);
  xbeam4 = beam.translateXyz(A, B, D);

  % Check all beams are equal
  target = xbeam1.getCoefficients;
  testCase.verifyThat(xbeam2.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'AB * beam does not match translateXyz');
  testCase.verifyThat(xbeam3.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'beam.translate(A, B) does not match translateXyz');
  testCase.verifyThat(xbeam4.getCoefficients, IsEqualTo(target, ...
    'Within', AbsoluteTolerance(tol)), ...
    'beam.translateXyz(A, B, D) does not match translateXyz');
end

function testUnevenTranslation(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  beam = ott.optics.vswf.bsc.PmGauss();
  beam.power = 1.0;
  dz = pi/2;
  tol = 1.0e-6;

  tbeam = beam.translateZ(dz, 'Nmax', beam.Nmax - 5);
  testCase.verifyThat(tbeam.Nmax, IsEqualTo(beam.Nmax - 5), ...
    'Translated beam does not have correct Nmax');

end

function testMakeBeamVectorEmpty(testCase)
% Check that make_beam_vector works for empty beams

  import matlab.unittest.constraints.IsEqualTo;

  a = [];
  b = [];
  nn = [];
  mm = [];
  
  [a1, b1] = ott.optics.vswf.bsc.Bsc.make_beam_vector(a, b, nn, mm);
  
  testCase.verifyThat([size(a1), size(b1)], IsEqualTo([0, 0, 0, 0]), ...
    'Size of beam vectors incorrect with implicit Nmax');
  
  [a2, b2] = ott.optics.vswf.bsc.Bsc.make_beam_vector(a, b, nn, mm, 1);
  
  testCase.verifyThat([size(a2), size(b2)], IsEqualTo([3, 0, 3, 0]), ...
    'Size of beam vectors incorrect with explicit Nmax');
end

function testMakeBeamVectorMulti(testCase)
% Check to make sure make_beam_vector functions with multiple beams
% with the same nn and mm indices.

  a = [1; 2; 3];
  b = [4; 5; 6];
  
  nn = [1; 2; 3];
  mm = [0; 0; 0];
  
  [a1, b1] = ott.optics.vswf.bsc.Bsc.make_beam_vector(a, b, nn, mm);
  [a2, b2] = ott.optics.vswf.bsc.Bsc.make_beam_vector(a+6, b+6, nn, mm);
  
  [ac, bc] = ott.optics.vswf.bsc.Bsc.make_beam_vector([a, a+6], [b, b+6], nn, mm);
  
  import matlab.unittest.constraints.IsEqualTo;
  
  testCase.verifyThat([ac(:, 1); bc(:, 1)], IsEqualTo([a1; b1]), ...
    'First beam doesn''t match');
  
  testCase.verifyThat([ac(:, 2); bc(:, 2)], IsEqualTo([a2; b2]), ...
    'Second beam doesn''t match');

end

function testSum(testCase)

  beam1 = ott.optics.vswf.bsc.PmGauss('polarisation', [0, 1i]);
  beam2 = ott.optics.vswf.bsc.PmGauss('polarisation', [1, 0]);
  beam3 = beam1 + beam2;
  
  beamU = beam1.append(beam2);
  testCase.verifyEqual(beamU.sum(), beam3, ...
    'beamU.sum() incorrect');
  testCase.verifyEqual(sum(beamU), beam3, ...
    'sum(beamU) incorrect');
  
  % Test array sum
  beamarr = [beam1, beam2];
  testCase.verifyEqual(sum(beamarr), beam3, ...
    'sum([beam1, beam2]) incorrect');
  
end

function testLargeTranslations(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  % For this translation the beam power should go to zero

  beam = ott.optics.vswf.bsc.PmGauss();
  beam.power = 1.0;
  tbeam = beam.translateXyz([300;0;0]);  % calls translate_z

  testCase.verifyThat(tbeam.power, IsEqualTo(0.0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam power does not drop to zero for large radial translations');

end
