function tests = tsetScattered
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

function testConstructTmatrix(testCase)

  ibeam = ott.beam.vswf.Bsc();
  tmatrix = ott.scat.vswf.Tmatrix();
  
  beam = ott.beam.vswf.Scattered.FromTmatrix(ibeam, tmatrix);
  
  testCase.verifyClass(beam, 'ott.optics.vswf.bsc.Scattered', 'cls');
  testCase.verifyEqual(beam.tmatrix, tmatrix, 'T');
  testCase.verifyEqual(beam.incident_beam, ibeam, 'ibeam');
  testCase.verifyEqual(beam.type, tmatrix.type, 'type');
  testCase.verifyEqual(beam.a, ibeam.a);
  testCase.verifyEqual(beam.b, ibeam.b);

end

function testTypeChangeWithIbeam(testCase)

  ibeam = ott.beam.vswf.Bsc();
  tmatrix = ott.scat.vswf.Tmatrix();
  
  beam = ott.beam.vswf.Scattered.FromTmatrix(ibeam, tmatrix);
  
  % TODO: Does this do what we want?
  
  beam.type = 'total';
  testCase.verifyEqual(beam.type, 'total', 'type T');
  
  beam.type = 'scattered';
  testCase.verifyEqual(beam.type, 'scattered', 'type S');
end

function testScatter(testCase)

  particle = ott.scat.vswf.Tmatrix.simple('sphere', 1.0, 'index_relative', 1.2);
  beam = ott.beam.vswf.PmGauss('power', 1.0);

  % Add a second beam
  beam0 = [beam, beam.translateZ(pi/2)];
  beam0.array_type = 'coherent';
  
  pos1 = [0;0;0];
  pos2 = [0;0;pi/2];
  
  sbeam0 = particle.scatter(beam0, 'position', [pos1, pos2]);
  
  beam1 = sum(beam0);
  beam1 = beam1.translateXyz(pos1);
  sbeam1 = particle.scatter(beam1);
  
  beam2 = sum(beam0);
  beam2 = beam2.translateXyz(pos2);
  sbeam2 = particle.scatter(beam2);
  
  testCase.verifyEqual(sbeam0, [sbeam1, sbeam2]);

end
