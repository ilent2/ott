function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

% Abstract class, no constructor

function testVisualise(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);

  h = figure();

  % Test visualise
  beam.visualise();

  % Test visualise far-field
  beam.visualiseFarfield();

  close(h);

end

function testCastArrayOfBeams(testCase)

  % Cast with two beams
  beam1 = ott.beam.abstract.Gaussian(1.0);
  beam2 = ott.beam.abstract.HermiteGaussian(1.0, 1, -1);
  beam_array = [beam1, beam2];

  beam = ott.beam.Beam(beam_array);
  testCase.assertClass(beam, 'ott.beam.Array');
  testCase.verifyEqual(beam.array_type, 'coherent');

end

function testCastArrayOfIncoherentArrays(testCase)

  % Cast with coherent and incoherent beams
  beam1 = ott.beam.abstract.Gaussian(1.0);
  beam2 = ott.beam.abstract.Incoherent([1, 1]);
  beam_array = [beam1, beam2];

  beam = ott.beam.Beam(beam_array);
  testCase.assertClass(beam, 'ott.beam.Array');
  testCase.verifyEqual(beam.array_type, 'array');

end

function testCastArray(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);

  beam = ott.beam.Array(beam);
  testCase.assertClass('ott.beam.Array');
  testCase.verifyEqual(beam.array_type, 'array');

end


function testCastCoherent(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);

  beam = ott.beam.Coherent(beam);
  testCase.assertClass('ott.beam.Array');
  testCase.verifyEqual(beam.array_type, 'coherent');

end

function testCastIncoherent(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);

  beam = ott.beam.Incoherent(beam);
  testCase.assertClass('ott.beam.Array');
  testCase.verifyEqual(beam.array_type, 'incoherent');

end

