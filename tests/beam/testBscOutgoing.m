function tests = testBscOutgoing
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testConstruct(testCase)

  index = 1.33;
  omega = 2*pi*3e8 ./ 532e-9;
  bsc = ott.bsc.Bsc([1;2;3], [4;5;6]);
  beam = ott.beam.BscOutgoing(bsc, 'index_medium', index, 'omega', omega);

  testCase.verifyEqual(beam.data, bsc, 'bsc');
  testCase.verifyEqual(beam.index_medium, index, 'index');
  testCase.verifyEqual(beam.omega, omega, 'omega');
end

function testFields(testCase)

  index = 1.33;
  omega = 2*pi*3e8 ./ 532e-9;
  bsc = ott.bsc.Bsc([1;2;3], [4;5;6]);
  beam = ott.beam.BscOutgoing(bsc, 'index_medium', index, 'omega', omega);

  rtp = [1;0;0];

  E = beam.efieldRtp(rtp);
  testCase.verifyInstanceOf(E, 'ott.utils.FieldVector', 'e');

  E = beam.hfieldRtp(rtp);
  testCase.verifyInstanceOf(E, 'ott.utils.FieldVector', 'h');
  
  % Test again with empty beam
  bsc = ott.bsc.Bsc(sparse(8, 1), sparse(8, 1));
  beam = ott.beam.BscOutgoing(bsc, 'index_medium', index, 'omega', omega);

  target = zeros(3, 1);
  E = beam.efieldRtp(rtp);
  testCase.verifyEqual(E.vxyz, target, 'e sp');
  
  E = beam.hfieldRtp(rtp);
  testCase.verifyEqual(E.vxyz, target, 'h sp');
end

