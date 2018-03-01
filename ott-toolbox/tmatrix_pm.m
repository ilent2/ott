function [T,T2] = tmatrix_pm(Nmax,k_medium,k_particle,shape,parameters)
% tmatrix_pm.m
% T-matrix, calculated by point-matching method
%
% Usage:
% T = tmatrix_pm(Nmax,k_medium,k_particle,shape,parameters);
% or
% [T,T2] = tmatrix_pm(Nmax,k_medium,k_particle,shape,parameters);
% where
% n_medium and n_particle are the wavenumbers in the
% surrounding medium and the scattering particle,
% shape and parameters describe the particle geometry
% (see shapesurface for more detail)
% T2 is the T matrix giving the internal field
%
% PACKAGE INFO

verbose = 0;

% Do we have a shape rotationally symmetric about the z-axis?
[dummy1,dummy2,rotational_symmetry] = shapesurface([],[],shape,parameters);
%rotational_symmetry = 0;

% If it's just a sphere, let's just find the Mie solution, which is, after all,
% analytical point-matching at a single point exploiting symmetry
%if rotational_symmetry == 2
%   [T,T2] = tmatrix_mie(Nmax,k_medium,k_particle,parameters(1));
%   return;
%end

total_orders = Nmax * (Nmax+2);

if rotational_symmetry
   ntheta = 4*(Nmax + 2);
   nphi = 1;
else
   ntheta = 1*(Nmax + 2);
   nphi = 2*(Nmax + 2)+1;
   ntheta = 2*(Nmax + 2);
   nphi = 3*(Nmax + 2)+1;
end

[theta,phi] = angulargrid(ntheta,nphi);

% Testing randomly distributed points
%xyz = randn(ceil(total_orders*50),3);
%[r,theta,phi] = xyz2rtp(xyz);

[r,normals] = shapesurface(theta,phi,shape,parameters);

npoints = length(theta);

% 3 vector components at each point, c/d,p/q coefficient per order
coeff_matrix = zeros(6*npoints,4*total_orders);
incident_wave_matrix = zeros(6*npoints,2*total_orders);

if rotational_symmetry
   T = sparse(2*total_orders,2*total_orders);
   T2 = sparse(2*total_orders,2*total_orders);
else
   T = zeros(2*total_orders,2*total_orders);
   T2 = zeros(2*total_orders,2*total_orders);
end

if verbose
   fprintf(1,'Generating data at points for each (n,m)\n');
end

for n = 1:Nmax
for m = -n:n

   % [M1,N1,M2,N2] = vswf(n,m,k_medium*r,theta,phi); % INCOMING-OUTGOING
   [M1,N1,dummy1,dummy2,M2,N2] = vswf(n,m,k_medium*r,theta,phi); % INCIDENT-SCATTERED
   [M3,N3] = vswf(n,m,k_particle*r,theta,phi,3);

   ci = combined_index(n,m);

   M1 = perpcomponent(M1,normals);
   N1 = perpcomponent(N1,normals);
   M2 = perpcomponent(M2,normals);
   N2 = perpcomponent(N2,normals);
   M3 = perpcomponent(M3,normals);
   N3 = perpcomponent(N3,normals);
   M1 = M1(:);
   N1 = N1(:);
   M2 = M2(:);
   N2 = N2(:);
   M3 = M3(:);
   N3 = N3(:);

   % 1 is outgoing field, 3 is particle field, 2 is incoming field
   coeff_matrix(:,ci) = - [ M1; N1 ];
   coeff_matrix(:,ci+total_orders) = - [ N1; M1 ];
   coeff_matrix(:,ci+2*total_orders) = [ M3; k_particle/k_medium*N3 ];
   coeff_matrix(:,ci+3*total_orders) = [ N3; k_particle/k_medium*M3 ];

   incident_wave_matrix(:,ci) = [ M2; N2 ];
   incident_wave_matrix(:,ci+total_orders) = [ N2; M2 ];

   if verbose
      fprintf(1,'.');
   end

end
end

if verbose
   fprintf(1,'done!\n');
   fprintf(1,'Calculating T-matrix columns\n');
end

for n = 1:Nmax
for m = -n:n

   ci = combined_index(n,m);

   if rotational_symmetry

      number_of_nm = 1 + Nmax - max(abs(m),1);
      nm_to_use = combined_index(max(abs(m),1):Nmax,ones(1,number_of_nm)*m);
      nm_to_use = [ nm_to_use nm_to_use+total_orders ];
      all_indices = [ nm_to_use  nm_to_use+2*total_orders ];

      incident_wave_vector = incident_wave_matrix(:,ci);
      Tcol = coeff_matrix(:,all_indices) \ incident_wave_vector;
      T(nm_to_use,ci) = Tcol(1:(2*number_of_nm),1);
      T2(nm_to_use,ci) = Tcol((1+2*number_of_nm):(4*number_of_nm),1);

      incident_wave_vector = incident_wave_matrix(:,ci+total_orders);
      Tcol = coeff_matrix(:,all_indices) \ incident_wave_vector;
      T(nm_to_use,ci+total_orders) = Tcol(1:(2*number_of_nm),1);
      T2(nm_to_use,ci+total_orders) = Tcol((1+2*number_of_nm):(4*number_of_nm),1);

   else

      incident_wave_vector = incident_wave_matrix(:,ci);
      Tcol = coeff_matrix \ incident_wave_vector;
      T(:,ci) = Tcol(1:2*total_orders,1);
      T2(:,ci) = Tcol((1+2*total_orders):4*total_orders,1);

      incident_wave_vector = incident_wave_matrix(:,ci+total_orders);
      Tcol = coeff_matrix \ incident_wave_vector;
      T(:,ci+total_orders) = Tcol(1:2*total_orders,1);
      T2(:,ci+total_orders) = Tcol((1+2*total_orders):4*total_orders,1);

   end

   if verbose
      fprintf(1,'.');
   end

end
end

if verbose
   fprintf(1,'done!\n');
end

return


