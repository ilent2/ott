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

  beam = ott.beam.abstract.Scattered(type, beam_data, ...
      'incident_beam', incident_beam);

  testCase.verifyEqual(beam.type, type, 'type');
  testCase.verifyEqual(beam.data_type, type, 'type');
  testCase.verifyEqual(beam.data, beam_data, 'beam_data');
  testCase.verifyEqual(beam.incident_beam, incident_beam, 'ibeam');

  % Verify dependent properties
  testCase.verifyEqual(beam.power, beam_data.power, 'power');
  testCase.verifyEqual(beam.medium, beam_data.medium, 'power');
  testCase.verifyEqual(beam.omega, beam_data.omega, 'power');

  % Verify total_beam and scattered_beam
  testCase.verifyEqual(beam.scattered_beam, beam_data, 'sbeam');
  testCase.verifyEqual(beam.total_beam, beam_data + incident_beam, 'tbeam');

end

function testContains(testCase)

  beam1 = ott.beam.abstract.Gaussian(1.0);
  array_beam = ott.beam.Array('incoherent', {beam1});
  beam = ott.beam.abstract.Scattered('scattered', array_beam);

  testCase.verifyEqual(beam.contains('incoherent'), true, '1incoherent');
  testCase.verifyEqual(beam.contains('coherent'), false, '1coherent');
  testCase.verifyEqual(beam.contains('array'), false, '1array');

end

function testCast(testCase)

  type = 'scattered';
  beam_data = ott.beam.abstract.Gaussian(1.0);
  incident_beam = ott.beam.abstract.Gaussian(2.0);

  abs_beam = ott.beam.abstract.Scattered(type, beam_data, ...
      'incident_beam', incident_beam);

  beam = ott.beam.Beam(abs_beam);
  testCase.verifyEqual(beam, beam_data);

end

