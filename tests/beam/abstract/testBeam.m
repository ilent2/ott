function tests = testBeam
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

% Abstract class, no constructor

function testContains(testCase)

  beam1 = ott.beam.abstract.Gaussian(1.0);
  beam2 = ott.beam.abstract.Gaussian(1.0);

  testCase.verifyEqual(beam1.arrayContains('incoherent'), false);
  testCase.verifyEqual(beam1.arrayContains('coherent'), false);
  testCase.verifyEqual(beam1.arrayContains('array'), false);

  beam = [beam1, beam2];
  testCase.verifyEqual(beam.arrayContains('incoherent'), false);
  testCase.verifyEqual(beam.arrayContains('coherent'), true);
  testCase.verifyEqual(beam.arrayContains('array'), false);

  beam = [beam1, beam2, ott.beam.Array('incoherent', [1, 1])];
  testCase.verifyEqual(beam.arrayContains('incoherent'), true);
  testCase.verifyEqual(beam.arrayContains('coherent'), true);
  testCase.verifyEqual(beam.arrayContains('array'), false);

end

function testArrayCast(testCase)

  beam1 = ott.beam.abstract.Gaussian(1.0);
  beam2 = ott.beam.abstract.Gaussian(1.0);

  cbeam = beam1.arrayCast(@ott.beam.paraxial.Gaussian);
  testCase.verifyClass(cbeam, 'ott.beam.paraxial.Gaussian');

  beam = [beam1, beam2];
  cbeam = beam.arrayCast(@ott.beam.paraxial.Gaussian);
  testCase.verifyClass(cbeam, 'ott.beam.Array');
  testCase.verifyEqual(cbeam.array_type, 'coherent');

  beam = [beam1, beam2, ott.beam.Array('incoherent', [1, 1])];
  testCase.verifyWarning(@() beam.arrayCast(@ott.beam.Beam), ...
    'ott:beam:abstract:Beam:not_coherent_array');

  S = warning('off', 'ott:beam:abstract:Beam:not_coherent_array');
  cbeam = beam.arrayCast(@ott.beam.Beam);
  warning(S);

  testCase.verifyClass(cbeam, 'ott.beam.Array');
  testCase.verifyEqual(cbeam.array_type, 'array');

end

function testVisualise(~)

  beam = ott.beam.abstract.Gaussian(1.0);

  h = figure();

  % Test visualise
  beam.visualise();

  % Test visualise far-field
  ott.beam.vswf.Bsc(beam, 'basis', 'incoming').visualiseFarfield();

  close(h);

end

function testOperators(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);

  testCase.verifyClass(-beam, 'ott.beam.abstract.Negative', 'uminus');

  testCase.verifyEqual(beam + beam, [beam, beam], 'plus');

  testCase.verifyEqual(beam - beam, [beam, -beam], 'minus');

end

function testCastArray(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);

  beam = ott.beam.Array(beam);
  testCase.assertClass(beam, 'ott.beam.Array');
  testCase.verifyEqual(beam.array_type, 'coherent');

end

