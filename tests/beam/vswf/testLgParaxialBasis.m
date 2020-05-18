function tests = testLgParaxialBasis
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../../');
end

% See testConsistency for default constructor

function testConstructOptional(testCase)

  lmodes = [0, 1, 2];
  pmodes = [0, 2, 4];
  polmodes = [-1, 1, 1];
  waist = 1.0;

  beam = ott.beam.vswf.LgParaxialBasis(waist, ...
    lmodes, pmodes, polmodes);

  testCase.verifyEqual(beam.waist, waist, 'waist');
  testCase.verifyEqual(beam.lmode, lmodes, 'lmodes');
  testCase.verifyEqual(beam.pmode, pmodes, 'amodes');
  testCase.verifyEqual(beam.polmode, polmodes, 'polmodes');
end

function testParaxialFields(testCase)

  beam = ott.beam.vswf.LgParaxialBasis();
  beam.waist = 1.0;
  
  theta = linspace(0, pi, 200).';
  lmode = 1;
  pmode = 5;
  polmode = 1;
  E = beam.paraxial_fields(theta, lmode, pmode, polmode);
  
  testCase.verifySize(E, [2*numel(theta), 1], 'sz');
  
%   figure();
%   plot(theta, abs(reshape(E, [], 2)));

end

