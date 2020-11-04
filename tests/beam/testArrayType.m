function tests = testArrayType
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)
  bsc = ott.bsc.Bsc(zeros(3, 2), ones(3, 2));
  beam = ott.beam.BscBeam(bsc, 'arrayType', 'coherent');
  testCase.verifyEqual(beam.arrayType, 'coherent');
end

function testCoherent(testCase)

  bsc = ott.bsc.Bsc(zeros(3, 2), ones(3, 2));
  beam = ott.beam.BscBeam(bsc, 'arrayType', 'coherent');
  
  im = beam.visNearfield();
  testCase.verifySize(im, [80, 80]);
end

function testIncoherent(testCase)

  bsc = ott.bsc.Bsc(zeros(3, 2), ones(3, 2));
  beam = ott.beam.BscBeam(bsc, 'arrayType', 'incoherent');
  
  im = beam.visFarfield();
  testCase.verifySize(im, [80, 80]);
end

function testArray(testCase)

  bsc = ott.bsc.Bsc(zeros(3, 2), ones(3, 2));
  beam = ott.beam.BscBeam(bsc, 'arrayType', 'array');
  
  im = beam.visFarfield();
  testCase.verifySize(im, [80, 80, 2]);
end
