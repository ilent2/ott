function tests = testIncecoefficients
  tests = functiontests(localfunctions);
end

function setupOnce(~)
  addpath('../../');
end

function testFromTheFile(testCase)

  %test in Bandres 2004, josaa-21-5-873.
  xi=3.;
  p=11;
  evalz=2;
  z=linspace(0,pi,100);

  [A_n, B_n, ~, ~] = ott.utils.incecoefficients(p, xi);

  C=zeros(size(A_n,1),100);
  
  delta=ceil(p/2)+1~=p/2+1;

  for ii=1:size(C,1)
      C(ii,:)=A_n(evalz,ii)*cos((2*(ii-1)+delta)*z);
  end

  h1 = figure();
  testCase.addTeardown(@() close(h1));
  plot(z./pi,sum(C(1:end,:),1),'b','linewidth',2);
  hold on
  plot([z(1),z(end)]/pi,[0,0],'k');
  hold off

  C=zeros(size(B_n,1),100);

  for ii=1:size(C,1)
      C(ii,:)=B_n(evalz,ii)*sin((2*(ii-1)+delta)*z);
  end

  h2 = figure();
  testCase.addTeardown(@() close(h2));
  plot(z./pi,sum(C(1:end,:),1),'b','linewidth',2);
  hold on
  plot([z(1),z(end)]/pi,[0,0],'k');
  hold off

end
