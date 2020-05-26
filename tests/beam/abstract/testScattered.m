function tests = testScattered
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructor(testCase)

  type = 'scattered';
  beam_data = ott.beam.abstract.Gaussian(1.0);
  incident_beam = ott.beam.abstract.Gaussian(2.0);
  particle = ott.scat.planewave.Plane([0;0;1], 1.5);

  beam = ott.beam.abstract.Scattered(type, beam_data, ...
      'incident_beam', incident_beam, 'particle', particle);

  testCase.verifyEqual(beam.type, type, 'type');
  testCase.verifyEqual(beam.beam_data, beam_data, 'beam_data');
  testCase.verifyEqual(beam.incident_beam, incident_beam, 'ibeam');
  testCase.verifyEqual(beam.particle, particle, 'particle');

  % Verify dependent properties
  testCase.verifyEqual(beam.power, beam_data.power, 'power');
  testCase.verifyEqual(beam.medium, beam_data.medium, 'power');
  testCase.verifyEqual(beam.omega, beam_data.omega, 'power');

  % TODO: Verify total_beam and scattered_beam

end

function testCast(testCase)

  type = 'scattered';
  beam_data = ott.beam.abstract.Gaussian(1.0);
  incident_beam = ott.beam.abstract.Gaussian(2.0);
  particle = ott.scat.planewave.Plane([0;0;1], 1.5);

  beam = ott.beam.abstract.Scattered(type, beam_data, ...
      'incident_beam', incident_beam, 'particle', particle);

  % Tests both Beam and Scattered
  beam = ott.beam.Beam(beam);

  testCase.assertClass(beam, 'ott.beam.Scattered');
  testCase.verifyEqual(beam.type, type, 'type');
  testCase.verifyEqual(beam.beam_data, beam_data, 'beam_data');
  testCase.verifyEqual(beam.incident_beam, incident_beam, 'ibeam');
  testCase.verifyEqual(beam.particle, particle, 'particle');

end

