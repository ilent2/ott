function [Ei_TE, Ei_TM] = E_inc_vswf(n,m,r_sp,k)
% Calculate the incident field for a single vswf mode
%
% [Ei_TE, Ei_TM] = E_inc_vswf(n,m,r_sp,k)
%
%    r_sp   Nx3 matrix of spherical coordinates [r, theta, phi]

N = size(r_sp, 1);
[M1,N1,M2,N2,M3,N3] = ott.utils.vswf(n,m,k*r_sp(:,1), r_sp(:,2), r_sp(:,3));
E_TE_sp = M3;
E_TM_sp = N3;

E_TE = zeros(3,N);
E_TM = zeros(3,N);
% convert to cartesian
% for j = 1:N
%   [E_TE(j,1),E_TE(j,2),E_TE(j,3)] = rtpv2xyzv(E_TE_sp(j,1),E_TE_sp(j,1),E_TE_sp(j,1),r_sp(j,1),r_sp(j,2),r_sp(j,3));
%   [E_TM(j,1),E_TM(j,2),E_TM(j,3)] = rtpv2xyzv(E_TM_sp(j,1),E_TM_sp(j,1),E_TM_sp(j,1),r_sp(j,1),r_sp(j,2),r_sp(j,3));
% end
for j = 1:N
  theta = r_sp(j,2);
  phi = r_sp(j,3);
  M = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
     cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
     -sin(phi) cos(phi) 0];
  E_TE(:,j) = M'*E_TE_sp(j,:)';
  E_TM(:,j) = M'*E_TM_sp(j,:)';
end


Ei_TE = zeros(3*N,1);
Ei_TM = zeros(3*N,1);

for j = 1:N % reformat into one column
  Ei_TE(3*(j-1) + 1) = E_TE(1,j);
  Ei_TE(3*(j-1) + 2) = E_TE(2,j);
  Ei_TE(3*(j-1) + 3) = E_TE(3,j);
  Ei_TM(3*(j-1) + 1) = E_TM(1,j);
  Ei_TM(3*(j-1) + 2) = E_TM(2,j);
  Ei_TM(3*(j-1) + 3) = E_TM(3,j);
end

