function tests = testBsc
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructDefault(testCase)

  beam = ott.beam.vswf.Bsc();
  
  testCase.verifyEqual(beam.a, [], 'a');
  testCase.verifyEqual(beam.b, [], 'b');
  testCase.verifyEqual(beam.basis, 'regular', 'basis');
  testCase.verifyEqual(beam.absdz, 0, 'dz');
  testCase.verifyEqual(beam.Nmax, 0, 'Nmax');
end

function testConstructAb(testCase)

  a = randn(8, 2);
  b = randn(8, 2);
  beam = ott.beam.vswf.Bsc(a, b);
  
  testCase.verifyEqual(beam.a, a, 'a');
  testCase.verifyEqual(beam.b, b, 'b');
  testCase.verifyEqual(beam.basis, 'regular', 'basis');
  testCase.verifyEqual(beam.absdz, 0, 'dz');
  testCase.verifyEqual(beam.Nmax, 2, 'Nmax');
  
  % Size related properties
  testCase.verifyEqual(numel(beam), 2, 'numel');
  testCase.verifySize(beam, [1, 2], 'sz');
  testCase.verifyNotEmpty(beam);
end

function testConstructBsc(testCase)

  a = randn(8, 2);
  b = randn(8, 2);
  beam0 = ott.beam.vswf.Bsc(a, b);
  beam1 = ott.beam.vswf.Bsc(beam0);
  
  testCase.verifyEqual(beam1, beam0);
end

function testTransformation(testCase)

  beam0 = ott.beam.vswf.Bsc(randn(8, 1), randn(8, 1));
  
  ab0 = beam0.getCoefficients();
  
  beam = beam0.translateXyz([0;0;1]);
  testCase.verifyEqual(beam.getCoefficients(), ab0, 'translate');
  testCase.verifyEqual(beam.position, [0;0;1], 'position');
  testCase.verifyClass(beam, 'ott.beam.vswf.Bsc');
  
  beamT = beam.applyTransformation();
  testCase.verifyEqual(beamT.position, [0;0;0], 'position');
  testCase.verifyClass(beamT, 'ott.beam.vswf.Bsc');
  
  beamS = beam0.applyTranslation([0;0;1]);
  testCase.verifyEqual(beamS.getCoefficients(), beamT.getCoefficients(), ...
    'AbsTol', 1.0e-15, 'cos');
  testCase.verifyEqual(beamS.position, [0;0;0], 'S position');
  testCase.verifyClass(beamS, 'ott.beam.vswf.Bsc');

end

function testTransformationArray(testCase)

  P = randn(3, 2);
  beam0 = ott.beam.vswf.Bsc(randn(8, 1), randn(8, 1));
  
  beam = beam0.translateXyz(P);
  testCase.verifySize(beam, [1, 2], 'translate sz');
  testCase.verifyClass(beam, 'ott.beam.vswf.Array', 'translate cls');
  
  beamT = beam.applyTransformation();
  testCase.verifySize(beamT, [1, 2], 'translate sz 2');
  testCase.verifyClass(beamT, 'ott.beam.vswf.Bsc', 'translate cls 2');
  
  beamS = beam0.applyTranslation(P);
  testCase.verifySize(beamS, [1, 2], 'Translation sz');
  testCase.verifyClass(beamS, 'ott.beam.vswf.Bsc', 'Translation cls');
  testCase.verifyEqual(beamS.getCoefficients(), beamT.getCoefficients(), ...
    'AbsTol', 1.0e-15, 'Array vs Direct');

end

function testTranslation(testCase)
  % Check that the translation functions run
  
  Nmax = 10;
  ab = randn(ott.utils.combined_index(Nmax, Nmax), 1);
  beam = ott.beam.vswf.Bsc(ab, ab);

  dz = pi/2;
  tol = 1.0e-6;

  tbeam1 = beam.applyZTranslation(dz);

  [~, AB] = beam.applyZTranslation(dz);
  tbeam2 = AB * beam;

  [~, A, B] = beam.applyZTranslation(dz);
  tbeam3 = [ A B ; B A ] * beam;
  tbeam4 = beam.applyTranslation('AB', {A; B});

  % Check all beams are equal
  target = tbeam1.getCoefficients;
  testCase.verifyEqual(tbeam2.getCoefficients, target, ...
    'AbsTol', tol, ...
    'AB * beam does not match translateZ');
  testCase.verifyEqual(tbeam3.getCoefficients, target, ...
    'AbsTol', tol, ...
    '[A, B; B, A] * beam does not match translateZ');
  testCase.verifyEqual(tbeam4.getCoefficients, target, ...
    'AbsTol', tol, ...
    'beam.translateZ(A, B, D) does not match translateZ');

  xbeam1 = beam.applyTranslation([dz; dz; 0]);

  [~, A, B, D] = beam.applyTranslation([dz; dz; 0]);
  xbeam4 = beam.applyTranslation('AB', {A; B; D});

  % Check all beams are equal
  target = xbeam1.getCoefficients;
  testCase.verifyEqual(xbeam4.getCoefficients, target, ...
    'AbsTol', tol, ...
    'beam.translateXyz(A, B, D) does not match translateXyz');
end

function testUnevenTranslation(testCase)

  Nmax = 20;
  ab = randn(ott.utils.combined_index(Nmax, Nmax), 1);
  beam = ott.beam.vswf.Bsc(ab, ab);
  
  tbeam = beam.applyZTranslation(1.0, 'Nmax', beam.Nmax - 5);
  testCase.verifyEqual(tbeam.Nmax, beam.Nmax - 5, ...
    'Translated beam does not have correct Nmax');

end

function testMakeBeamVectorEmpty(testCase)

  a = [];
  b = [];
  nn = [];
  mm = [];
  
  bsc = ott.beam.vswf.Bsc.FromDenseBeamVectors(a, b, nn, mm);
  testCase.verifyEmpty(bsc);
end

function testMakeBeamVectorMulti(testCase)
% Check to make sure make_beam_vector functions with multiple beams
% with the same nn and mm indices.

  a = [1; 2; 3];
  b = [4; 5; 6];
  
  nn = [1; 2; 3];
  mm = [0; 0; 0];
  
  bsc1 = ott.beam.vswf.Bsc.FromDenseBeamVectors(a, b, nn, mm);
  bsc2 = ott.beam.vswf.Bsc.FromDenseBeamVectors(a+6, b+6, nn, mm);
  
  bsc3 = ott.beam.vswf.Bsc.FromDenseBeamVectors([a, a+6], [b, b+6], nn, mm);
  
  testCase.verifyEqual(bsc3(1), bsc1, ...
    'First beam doesn''t match');
  testCase.verifyEqual(bsc3(2), bsc2, ...
    'Second beam doesn''t match');
end

function testSum(testCase)

  Nmax = 5;
  ab1 = randn(ott.utils.combined_index(Nmax, Nmax), 1);
  ab2 = randn(ott.utils.combined_index(Nmax, Nmax), 1);
  beam1 = ott.beam.vswf.Bsc(ab1, ab1);
  beam2 = ott.beam.vswf.Bsc(ab2, ab2);
  
  beam3 = beam1 + beam2;
  testCase.assertSize(beam3, [1, 1], 'sz');
  
  beamU = [beam1, beam2];
  testCase.verifyEqual(beamU.sum(), beam3, ...
    'beamU.sum() incorrect');
  testCase.verifyEqual(sum(beamU), beam3, ...
    'sum(beamU) incorrect');
  
end

function testLargeTranslations(testCase)
  % For this translation the beam power should go to zero

  beam = ott.beam.vswf.Gaussian();
  beam.power = 1.0;
  tbeam = beam.applyTranslation([6000;0;0]);  % calls translate_z

  testCase.verifyEqual(tbeam.power, 0.0, ...
      'AbsTol', 1.0e-6, ...
      'Beam power does not drop to zero for large radial translations');

end
