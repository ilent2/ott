function tests = testLgParaxial
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testSingleMode(testCase)

  lmode = 1;
  pmode = 2;
  polmode = 1;
  waist = 1;
  beam = ott.bsc.LgParaxial(waist, lmode, pmode, polmode);

  testCase.verifyEqual(beam.lmode, lmode, 'l');
  testCase.verifyEqual(beam.pmode, pmode, 'p');
  testCase.verifyEqual(beam.polmode, polmode, 'pol');
  testCase.verifyEqual(beam.waist, waist, 'w');

end

function testMultipleModes(testCase)

  lmode = [-1, 0, 1];
  pmode = 2;
  polmode = 1;
  waist = 1;
  beam = ott.bsc.LgParaxial(waist, lmode, pmode, polmode);

  testCase.verifyEqual([beam.lmode], lmode, 'l');
  testCase.verifyEqual([beam.pmode], [1,1,1]*pmode, 'p');
  testCase.verifyEqual([beam.polmode], [1,1,1]*polmode, 'pol');
  testCase.verifyEqual([beam.waist], [1,1,1]*waist, 'w');

end

function testFromIg(testCase)

  waist = 1;
  lmode = 1;
  porder = 1;
  parity = 'even';
  ellipticity = 1;
  polbasis = 'cartesian';
  polfield = [1, 1i];
  bsc = ott.bsc.LgParaxial.FromIgMode(waist, lmode, porder, ...
    parity, ellipticity, polbasis, polfield);
  
  testCase.verifyClass(bsc, 'ott.bsc.Bsc');

end

function testFromHg(testCase)

  waist = 1;
  mmode = 1;
  nmode = 1;
  polbasis = 'cartesian';
  polfield = [1, 1i];
  
  bsc = ott.bsc.LgParaxial.FromHgMode(waist, mmode, nmode, ...
      polbasis, polfield);
  
  testCase.verifyClass(bsc, 'ott.bsc.Bsc');

end

function testFromLg(testCase)

  waist = 1;
  lmode = 1;
  pmode = 1;
  polbasis = 'cartesian';
  polfield = [1, 1i];
  
  bsc = ott.bsc.LgParaxial.FromLgMode(waist, lmode, pmode, ...
      polbasis, polfield);
  
  testCase.verifyClass(bsc, 'ott.bsc.Bsc');

end
