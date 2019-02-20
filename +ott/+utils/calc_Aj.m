function Aj = calc_Aj(k,r,alph,j,blockdiag)
% calculates a 3 X 3N block comprising N number of 3 X 3 Green's tensors  

if nargin < 5
  blockdiag = 1;
end

%   rk_to_rj = r(jj,:)-r(kk,:);
%   rjk = norm(rk_to_rj);
%   rjk_x = rk_to_rj(1)/rjk;
%   rjk_y = rk_to_rj(2)/rjk;
%   rjk_z = rk_to_rj(3)/rjk;
  [N,col] = size(r);
  rk_to_rj = repmat(r(j,:),N,1)-r;
  rjk=sqrt(sum(rk_to_rj.^2,2));
  rjk_x = rk_to_rj(:,1)./rjk;                       rjk_x(j) = 0;
  rjk_y = rk_to_rj(:,2)./rjk;                       rjk_y(j) = 0;
  rjk_z = rk_to_rj(:,3)./rjk;                       rjk_z(j) = 0;
%   beta = (1 - (k*rjk)^-2 + i*(k*rjk)^-1);
%   gamma = -(1 - 3*(k*rjk)^-2 + i*3*(k*rjk)^-1);
  B = (1 - (k*rjk).^-2 + i*(k*rjk).^-1);            B(j) = 0;
  G = -(1 - 3*(k*rjk).^-2 + i*3*(k*rjk).^-1);       G(j) = 0;
  Aj = zeros(3*N,3); % vertical at first

%   Ajk = -exp(i*k*rjk)/rjk*k^2*...
%     [(beta + gamma*rjk_x^2) (gamma*rjk_x*rjk_y)    (gamma*rjk_x*rjk_z);
%     (gamma*rjk_y*rjk_x)    (beta + gamma*rjk_y^2) (gamma*rjk_y*rjk_z);
%     (gamma*rjk_z*rjk_x)    (gamma*rjk_z*rjk_y)    (beta + gamma*rjk_z^2)];
  C = -exp(i*k*rjk)./rjk*k^2;                       C(j) = 0;
  Aj([1:3:3*N],1) = (C.*(B + G.*rjk_x.^2));
  Aj([2:3:3*N],1) = (C.*(G.*rjk_x.*rjk_y));
  Aj([3:3:3*N],1) = (C.*(G.*rjk_x.*rjk_z));
  Aj([1:3:3*N],2) = (C.*(G.*rjk_y.*rjk_x));
  Aj([2:3:3*N],2) = (C.*(B + G.*rjk_y.^2));
  Aj([3:3:3*N],2) = (C.*(G.*rjk_y.*rjk_z));
  Aj([1:3:3*N],3) = (C.*(G.*rjk_z.*rjk_x));
  Aj([2:3:3*N],3) = (C.*(G.*rjk_z.*rjk_y));
  Aj([3:3:3*N],3) = (C.*(B + G.*rjk_z.^2));
  Aj = conj(Aj'); % now tranpose and unconjugate it

  % block diagonal - inverse polarizability tensor
  if blockdiag
    Aj(1,(j-1)*3+1) = 1/alph((j-1)*3+1);
    Aj(2,(j-1)*3+2) = 1/alph((j-1)*3+2);
    Aj(3,(j-1)*3+3) = 1/alph((j-1)*3+3);
  end
  
