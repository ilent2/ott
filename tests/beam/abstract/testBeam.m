function tests = testBessel
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testVisualise(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);

  h = figure();

  % Test visualise
  beam.visualise();

  % Test visualise far-field
  beam.visualiseFarfield();

  close(h);

end

function testCastBeam(testCase)

  beam1 = ott.beam.abstract.Gaussian(1.0);
  beam2 = ott.beam.abstract.HermiteGaussian(1.0, 1, -1);
  beam_array = [beam1, beam2];

  beam = ott.beam.Beam(beam_array);
  testCase.verifyClass('ott.beam.Coherent');

end

function testCasts(testCase)

  beam = ott.beam.abstract.Gaussian(1.0);

  beam = ott.beam.Array(beam);
  testCase.verifyClass('ott.beam.Array');

  beam = ott.beam.Coherent(beam);
  testCase.verifyClass('ott.beam.Coherent');

  beam = ott.beam.Incoherent(beam);
  testCase.verifyClass('ott.beam.Incoherent');

end

