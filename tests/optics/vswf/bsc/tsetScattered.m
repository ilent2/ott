function tests = tsetScattered
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../../');
end

function testConstructAbType(testCase)

  type = 'scattered';
  a = [];
  b = [];
  
  beam = ott.optics.vswf.bsc.Scattered(a, b, type);
  
  testCase.verifyClass(beam, 'ott.optics.vswf.bsc.Scattered', 'cls');
  testCase.verifyEqual(beam.tmatrix, [], 'T');
  testCase.verifyEqual(beam.incident_beam, [], 'ibeam');
  testCase.verifyEqual(beam.type, type, 'type');
  testCase.verifyEqual(beam.a, a, 'a');
  testCase.verifyEqual(beam.b, b, 'b');
  
end

function testConstructBscType(testCase)

  bsc = ott.optics.vswf.bsc.Bsc();
  type = 'total';
  
  beam = ott.optics.vswf.bsc.Scattered(bsc, type);
  
  testCase.verifyClass(beam, 'ott.optics.vswf.bsc.Scattered', 'cls');
  testCase.verifyEqual(beam.tmatrix, [], 'T');
  testCase.verifyEqual(beam.incident_beam, [], 'ibeam');
  testCase.verifyEqual(beam.type, type, 'type');
  testCase.verifyEqual(beam.a, bsc.a, 'a');
  testCase.verifyEqual(beam.b, bsc.b, 'b');
end

function testConstructTmatrix(testCase)

  ibeam = ott.optics.vswf.bsc.Bsc();
  tmatrix = ott.optics.vswf.tmatrix.Tmatrix();
  
  beam = ott.optics.vswf.bsc.Scattered(ibeam, tmatrix);
  
  testCase.verifyClass(beam, 'ott.optics.vswf.bsc.Scattered', 'cls');
  testCase.verifyEqual(beam.tmatrix, tmatrix, 'T');
  testCase.verifyEqual(beam.incident_beam, ibeam, 'ibeam');
  testCase.verifyEqual(beam.type, tmatrix.type, 'type');
  testCase.verifyEqual(beam.a, ibeam.a);
  testCase.verifyEqual(beam.b, ibeam.b);

end

function testTypeChangeWithIbeam(testCase)

  ibeam = ott.optics.vswf.bsc.Bsc();
  tmatrix = ott.optics.vswf.tmatrix.Tmatrix();
  
  beam = ott.optics.vswf.bsc.Scattered(ibeam, tmatrix);
  
  beam.type = 'total';
  testCase.verifyEqual(beam.type, 'total', 'type T');
  
  beam.type = 'scattered';
  testCase.verifyEqual(beam.type, 'scattered', 'type S');
end
