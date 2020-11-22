function tests = testArrays
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
  
  % Generate beams for testing
  testCase.TestData.beam1 = ott.beam.Gaussian();
end

function testConstruct(testCase)

  beam1 = ott.beam.BscBeam();
  beam2 = ott.beam.BscBeam();
  beams = [beam1, beam2];
  
  beamA = ott.beam.Array(beams, 'arrayType', 'array');
  
  testCase.verifyEqual(beamA.data, beams, 'data');
  testCase.verifyEqual(beamA.arrayType, 'array', 'type');
end

function testEmptyArray(testCase)

  beam = ott.beam.Beam.empty(1, 0);
  beamA = ott.beam.Array(beam, 'arrayType', 'array');
  
  testCase.verifyEqual(beamA.index_medium, nan, 'index');
  testCase.verifyEqual(beamA.omega, nan, 'omega');
  
end

function testArrayOfArrays(testCase)

  beam = ott.beam.Empty('index_medium', 1.33);
  beamA = ott.beam.Array([beam, beam], 'arrayType', 'coherent');
  beamB = ott.beam.Array([beamA, beam], 'arrayType', 'array');
  
  testCase.verifyEqual(beamB.index_medium, beam.index_medium, 'index');

end

function testCoherent(testCase)

  beam(1) = testCase.TestData.beam1;
  beam(2) = testCase.TestData.beam1;
  beam(1).position = [-1;0;0]*1e-6;
  beam(2).position = [1;0;0]*1e-6;
  sz = [10,10];

  e1 = beam(1).visNearfield('range', [1,1]*5e-6, ...
    'field', 'Ex', 'size', sz);
  e2 = beam(2).visNearfield('range', [1,1]*5e-6, ...
    'field', 'Ex', 'size', sz);
  es = e1 + e2;

  cbeam = ott.beam.Array(beam, 'arrayType', 'coherent');
  e = cbeam.visNearfield('range', [1, 1]*5e-6, 'field', 'Ex', 'size', sz);

  testCase.verifyEqual(e, es, 'Ex');

end

function testIncoherent(testCase)

  beam(1) = testCase.TestData.beam1;
  beam(2) = testCase.TestData.beam1;
  beam(1).position = [-1;0;0]*1e-6;
  beam(2).position = [1;0;0]*1e-6;
  sz = [10,10];

  e1 = beam(1).visNearfield('range', [1,1]*5e-6, ...
    'field', 'irradiance', 'size', sz);
  e2 = beam(2).visNearfield('range', [1,1]*5e-6, ...
    'field', 'irradiance', 'size', sz);
  es = e1 + e2;

  cbeam = ott.beam.Array(beam, 'arrayType', 'incoherent');
  e = cbeam.visNearfield('range', [1, 1]*5e-6, ...
    'field', 'irradiance', 'size', sz);

  testCase.verifyEqual(e, es, 'irradiance');

end

function testPosition(testCase)

  offset = [0;1;0]*1e-6;

  beam(1) = testCase.TestData.beam1;
  beam(2) = testCase.TestData.beam1;
  beam(1).position = [-1;0;0]*1e-6;
  beam(2).position = [1;0;0]*1e-6;
  sz = [10,10];
  
  cbeam1 = ott.beam.Array(beam, 'arrayType', 'incoherent');
  cbeam1.position = offset;

  beam = beam.translateXyz(offset);
  cbeam2 = ott.beam.Array(beam, 'arrayType', 'incoherent');

  e1 = cbeam1.visNearfield('range', [1,1]*1e-6, 'size', sz);
  e2 = cbeam2.visNearfield('range', [1,1]*1e-6, 'size', sz);

  testCase.verifyEqual(e1, e2, 'position incorrect');

end

function testFieldCoverage(testCase)

  beam = testCase.TestData.beam1;
  cbeam = ott.beam.Array(beam, 'arrayType', 'incoherent');
  
  cbeam.efarfield([0;0;0]);
  cbeam.hfarfield([0;0;0]);
  cbeam.ehfarfield([0;0;0]);
  cbeam.eparaxial([0;0]);
  cbeam.hparaxial([0;0]);
  cbeam.ehparaxial([0;0]);

end
