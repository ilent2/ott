function tests = testVswf
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testMultipleOutputs(testCase)

  npts = 5;
  n = 3;
  m = -n:n;
  k_medium = 1.33;
  r = rand(npts, 1);
  theta = pi*rand(npts, 1);
  phi = 2*pi*rand(npts, 1);
  [M1,N1,M2,N2,M3,N3] = ott.utils.vswf(n,m,k_medium*r,theta,phi);
  
  testCase.assertSize(M1, [5, 3*numel(m)], 'sz M1');
  testCase.assertSize(N1, [5, 3*numel(m)], 'sz N1');
  testCase.assertSize(M2, [5, 3*numel(m)], 'sz M2');
  testCase.assertSize(N2, [5, 3*numel(m)], 'sz N2');
  testCase.assertSize(M3, [5, 3*numel(m)], 'sz M3');
  testCase.assertSize(N3, [5, 3*numel(m)], 'sz N3');
  
  for ii = 1:numel(m)
    [M1t,N1t,M2t,N2t,M3t,N3t] = ott.utils.vswf(n,m(ii),k_medium*r,theta,phi);
    idx = (0:2)*numel(m) + ii;
    
    testCase.verifyEqual(M1(:, idx), M1t, 'mt1');
    testCase.verifyEqual(N1(:, idx), N1t, 'nt1');
    testCase.verifyEqual(M2(:, idx), M2t, 'mt2');
    testCase.verifyEqual(N2(:, idx), N2t, 'nt2');
    testCase.verifyEqual(M3(:, idx), M3t, 'mt3');
    testCase.verifyEqual(N3(:, idx), N3t, 'nt3');
  end

end
