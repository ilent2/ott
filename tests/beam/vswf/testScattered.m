function tests = testScattered()
  % Test scattered beam
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

function testConstructAbType(testCase)

  type = 'scattered';
  a = [];
  b = [];

  beam = ott.beam.vswf.Scattered(a, b, 'type', type);
  testCase.verifyClass(beam, 'ott.beam.vswf.Scattered', 'cls');
  testCase.verifyEmpty(beam.tmatrix, 'T');
  testCase.verifyEmpty(beam.incident_beam, 'ibeam');
  testCase.verifyEqual(beam.type, type, 'type');
  testCase.verifyEqual(beam.a, a, 'a');
  testCase.verifyEqual(beam.b, b, 'b');
end

function testConstructBscType(testCase)

  bsc = ott.beam.vswf.Bsc();
  type = 'total';
  beam = ott.beam.vswf.Scattered(bsc, 'type', type);

  testCase.verifyClass(beam, 'ott.beam.vswf.Scattered', 'cls');
  testCase.verifyEmpty(beam.tmatrix, 'T');
  testCase.verifyEmpty(beam.incident_beam, 'ibeam');
  testCase.verifyEqual(beam.type, type, 'type');
  testCase.verifyEqual(beam.a, bsc.a, 'a');
  testCase.verifyEqual(beam.b, bsc.b, 'b');
end

