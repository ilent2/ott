function tests = testArrayAdapter
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../../');
end

function testConstruct(testCase)

  array = ott.beam.Array.empty('array_type', 'array');
  beam = ott.beam.abstract.ArrayAdapter(array);

  testCase.verifyEqual(beam.data, array);

end

