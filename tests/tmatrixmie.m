function tests = tmatrixmie
  tests = functiontests(localfunctions);
end

function testConstruct(testCase)

  addpath('../');

  % Simple sphere
  Tsimple = ott.TmatrixMie(1.0, 'index_relative', 1.2);

  % Layered sphere
  Tsimple = ott.TmatrixMie([0.5, 1.0], 'index_relative', [1.4, 1.2]);

end

