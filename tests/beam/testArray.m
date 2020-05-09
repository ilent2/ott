function tests = testArray
  % Test for beam combination functionality
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testPlusSameArrayType(testCase)

  beam1 = ott.beam.PlaneWave('direction', randn(3, 3));
  beam2 = ott.beam.PlaneWave('direction', randn(3, 4));
  
  testCase.verifyEqual(beam1.array_type, 'coherent', 'default value');
  testCase.verifyEqual(size(beam1), [1, 3], 'original size');
  
  % Adding a beam to a coherent array should just extend the array
  new_beam = beam1 + beam2;
  testCase.verifyEqual(size(new_beam), [1, 7], 'size');
  testCase.verifyEqual(new_beam.array_type, 'coherent', 'new type');
  
  % Adding beams to an array type should create a coherent type
  beam1.array_type = 'array';
  new_beam = beam1 + beam2;
  testCase.verifyEqual(new_beam.array_type, 'coherent', 'array new type');
  testCase.verifyEqual(size(new_beam), [1, 2], 'array size');
  testCase.verifyEqual(size(new_beam(1)), [1, 3], 'sub-array size');
  
  % Adding beams to a incoherent array should cause an error
  beam1.array_type = 'incoherent';
  testCase.verifyError(@() beam1 + beam2, ...
    'ott:beam:utils:ArrayType:plus_incoherent');

end

function testDifferentType(testCase)

  beam1 = ott.beam.PlaneWave();
  beam2 = ott.beam.paraxial.Gaussian(1.0);
  
  new_beam = beam1 + beam2;
  testCase.verifyEqual(size(new_beam), [1, 2], 'different types');
  testCase.verifyEqual(new_beam.array_type, 'coherent', 'array new type');

end

function testArrayCat(testCase)

  beam1 = ott.beam.PlaneWave();
  beam2 = ott.beam.paraxial.Gaussian(1.0);
  
  new_beam = [beam1, beam2; beam2, beam1];
  testCase.verifyEqual(size(new_beam), [2, 2], 'square :)');
  
end

function testArrayCatIncoherentError(testCase)

  beam1 = ott.beam.PlaneWave();
  beam2 = ott.beam.paraxial.Gaussian(1.0);
  beam = [beam1, beam2];
  
  beam.array_type = 'incoherent';
  
  testCase.verifyError(@() ott.beam.utils.ArrayType.AutoArray('coherent', [1, 2], beam, beam), ...
    'ott:beam:utils:ArrayType:coherent_with_incoherent');

end

function testLogicalSelection(testCase)

  beam1 = ott.beam.PlaneWave();
  beam2 = ott.beam.paraxial.Gaussian(1.0);
  beam = [beam1, beam2];
  
  target1 = ott.beam.Array('array', [1, 1], beam2);
  target3 = ott.beam.Array('array', [0, 0]);

  testCase.verifyEqual(beam([false, true]), target1, 'beam2');
  testCase.verifyEqual(beam([true, true]), beam, 'beam');
  testCase.verifyEqual(beam([false, false]), target3, 'empty');

end

function testLogicalSubSelection(testCase)

  plane = ott.beam.PlaneWave('origin', randn(3, 5));
  beam = ott.beam.Array('array', [1, 2], plane, plane);
  larray = [false, true, false, false, true];
  farray = false(1, 5);
  
  testCase.verifyEqual(beam(1), plane, '()');
  testCase.verifyEqual(beam.beams{1}, plane, '.beams{}');
  testCase.verifyEqual(beam.beams{1}(larray), plane(larray), 'larray');
  testCase.verifyEqual(beam.beams{1}(farray), plane(farray), 'farray');

end

function testEmptyCat(testCase)

  empty = ott.beam.abstract.Empty();
  plane = ott.beam.PlaneWave('origin', randn(3, 5));
  
  beam = [empty, plane];
  testCase.verifyEqual(beam, plane, 'empty cat');
  
  array = ott.beam.Array('array', [1, 3]);
  testCase.verifyEqual(array(1), empty, 'empty array');
  testCase.verifyEqual([array(1), plane], plane, 'more empty cat');
end

