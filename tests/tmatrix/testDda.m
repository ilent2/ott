function tests = testDda
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testSphereFromShape(testCase)

  radius = 0.2;
  nrel = 1.2;
  shape = ott.shape.Sphere(radius);
  Tdda = ott.tmatrix.Dda.FromShape(shape, ...
      'relative_index', nrel, 'spacing', 1/41);
  Tmie = ott.tmatrix.Mie.FromShape(shape, 'relative_index', nrel);

  testCase.assertEqual(Tdda.Nmax, Tmie.Nmax, 'Nmax');
  testCase.verifyEqual(Tdda.data, Tmie.data, 'AbsTol', 1e-5, 'data');

end

function testDipoleSphere(testCase)

  xyz = [0;0;0];
  nrel = 1.2;
  
  radius = 0.01;
  d = (4*pi/3).^(1/3) .* radius;
  
  dda = ott.tmatrix.dda.Dda(xyz, ...
    ott.tmatrix.dda.polarizability.CM(d, nrel)*eye(3));
  Tdda = ott.tmatrix.Dda(dda, 'Nmax', 1);
  
  shape = ott.shape.Sphere(radius);
  Tmie = ott.tmatrix.Mie.FromShape(shape, ...
      'relative_index', nrel, 'Nmax', 1);
  
  Tdda_diag = diag(Tdda);
  Tmie_diag = diag(Tmie);
  testCase.verifyEqual(Tdda_diag, Tmie_diag, 'diag');
%   testCase.verifyEqual(Tdda_diag(4:end), Tmie_diag(4:end), ...
%     'RelTol', 0.001, 'Upper Mie coefficeints not equal');
%   testCase.verifyEqual(Tdda.data, full(Tmie.data), ...
%     'AbsTol', 1e-6, 'Upper Mie coefficeints not equal');

end

function testSphereNearfield(testCase)
  % Check that near-field matches far-field

  radius = 0.2;
  nrel = 1.2;
  shape = ott.shape.Sphere(radius);
  Nmax = ott.utils.ka2nmax(2*pi*radius);
  pmrtp = ott.tmatrix.Dda.DefaultPmrtp(Nmax, 'radius', 8, ...
    'xySymmetry', true, 'zRotSymmetry', 0);
  
  TddaNear = ott.tmatrix.Dda.FromShape(shape, ...
      'relative_index', nrel, 'spacing', 1/35, ...
      'pmrtp', pmrtp, 'ci', 1);
    
  TddaFar = ott.tmatrix.Dda.FromShape(shape, ...
      'relative_index', nrel, 'spacing', 1/35, 'ci', 1);
  
  testCase.verifyTrue(all(all(TddaNear.data(:, [2,3,5,6]) == 0)), 'zero columns');
  testCase.verifyEqual(TddaNear.data, TddaFar.data, 'AbsTol', 1e-6, 'data');

end
