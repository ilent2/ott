function [T,T2] = tmatrix_pm_cube(Nmax,Nmax_medium,Nmax_particle,k_medium,k_particle,sidelength)
% tmatrix_pm_cube.m
% T-matrix, calculated by point-matching method
%
% Usage:
% T =
% tmatrix_pm_cube(Nmax,Nmax_medium,Nmax_particle,k_medium,k_particle,sidelength)
% or
% [T,T2] =
% tmatrix_pm_cube(Nmax,Nmax_medium,Nmax_particle,k_medium,k_particle,sidelength)
% where
% n_medium and n_particle are the wavenumbers in the
% surrounding medium and the scattering particle,
% shape and parameters describe the particle geometry
% (see shapesurface for more detail)
% T2 is the T matrix giving the internal field
%
% PACKAGE INFO

verbose = 0;
shape=4;

% Do we have a shape rotationally symmetric about the z-axis?
[dummy1,dummy2,rotational_symmetry] = shapesurface([],[],shape,sidelength);
%rotational_symmetry = 0;

% If it's just a sphere, let's just find the Mie solution, which is, after all,
% analytical point-matching at a single point exploiting symmetry
%if rotational_symmetry == 2
%   [T,T2] = tmatrix_mie(Nmax,k_medium,k_particle,parameters(1));
%   return;
%end

% calculates the size of the matrices
total_orders = Nmax * (Nmax+2);
total_orders_medium = Nmax_medium * (Nmax_medium+2);
total_orders_particle = Nmax_particle * (Nmax_particle+2);

% calculates the number of points for the grid
ntheta = 2*(Nmax_particle + 2)+1;
nphi = 3*(Nmax_particle + 2)+1;
   
[theta,phi] = angulargrid(round(ntheta),round(nphi));

% Testing randomly distributed points
%xyz = randn(ceil(total_orders*50),3);
%[r,theta,phi] = xyz2rtp(xyz);
[r,normals] = shapesurface(theta,phi,shape,sidelength);
npoints = length(theta);

% 3 vector components at each point, c/d,p/q coefficient per order
coeff_matrix = zeros(6*npoints,2*total_orders_particle+2*total_orders_medium);
incident_wave_matrix = zeros(6*npoints,2*total_orders_medium);

% Sets up the T matrices to be sparse
T = sparse(2*total_orders,2*total_orders);
T2 = sparse(2*total_orders,2*total_orders);


if verbose
   fprintf(1,'Generating data at points for each (n,m)\n');
end

if verbose
    h = waitbar(0,'Generating data at points for each (n,m)');
end

for n = 1:Nmax_medium
for m = -n:n

   % [M1,N1,M2,N2] = vswf(n,m,k_medium*r,theta,phi); % INCOMING-OUTGOING
   [M1,N1,dummy1,dummy2,M2,N2] = vswf(n,m,k_medium*r,theta,phi); % INCIDENT-SCATTERED
   [M3,N3] = vswf(n,m,k_particle*r,theta,phi,3);

   ci = combined_index(n,m);
   % This bit is just for the waitbar
   ci_max = combined_index(Nmax_medium,Nmax_medium);

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
   coeff_matrix(:,ci+total_orders_medium) = - [ N1; M1 ];
   coeff_matrix(:,ci+2*total_orders_medium) = [ M3; k_particle/k_medium*N3 ];
   coeff_matrix(:,ci+(2*total_orders_medium+total_orders_particle)) = [ N3; k_particle/k_medium*M3 ];

   size(coeff_matrix);
   
   incident_wave_matrix(:,ci) = [ M2; N2 ];
   incident_wave_matrix(:,ci+total_orders_medium) = [ N2; M2 ];
   
   size(incident_wave_matrix);

   if verbose
      fprintf(1,'.');
   end

   if verbose
       waitbar(ci/ci_max,h)
   end
end
end
if verbose
    close(h)
end

if verbose
    h = waitbar(0,'Calculating T-matrix columns...');
end
if verbose
   fprintf(1,'done!\n');
   fprintf(1,'Calculating T-matrix columns\n');
end

for n = 1:Nmax_medium
for m = -n:n

   ci = combined_index(n,m);    
   % This bit is just for the waitbar
   ci_max = combined_index(Nmax_medium,Nmax_medium);

      [n_medium_same m_medium_same ci_medium_same] = nm_cube(n,m,Nmax_medium);
      [n_cube_same m_cube_same ci_cube_same] = nm_cube(n,m,Nmax_particle);
      [n_medium_opp m_medium_opp ci_medium_opp] = nm_cube(n+1,m,Nmax_medium);
      [n_cube_opp m_cube_opp ci_cube_opp] = nm_cube(n+1,m,Nmax_particle);
      
      % TE
      indices_medium = [ ci_medium_same' (ci_medium_opp+total_orders_medium)' ];
      indices_cube = [ (ci_cube_same+2*total_orders_medium)' (ci_cube_opp+2*total_orders_medium+total_orders_particle)'];
      all_indices = [ indices_medium indices_cube ];
      out_medium = [ ci_medium_same' (ci_medium_opp+total_orders)' ];
      out_cube = [ ci_cube_same' (ci_cube_opp+total_orders)'];
      incident_wave_vector = incident_wave_matrix(:,ci);
      Tcol = coeff_matrix(:,all_indices) \ incident_wave_vector;
      T(out_medium,ci) = Tcol(1:length(indices_medium),1);
      T2(out_cube,ci) = Tcol((1+length(indices_medium)):length(all_indices),1);

      % TM
      indices_medium = [ ci_medium_opp' (ci_medium_same+total_orders_medium)' ];
      indices_cube = [ (ci_cube_opp+2*total_orders_medium)' (ci_cube_same+2*total_orders_medium+total_orders_particle)'];
      all_indices = [ indices_medium indices_cube ];
      out_medium = [ ci_medium_opp' (ci_medium_same+total_orders)' ];
      out_cube = [ ci_cube_opp' (ci_cube_same+total_orders)'];
      incident_wave_vector = incident_wave_matrix(:,ci+total_orders_medium);
      Tcol = coeff_matrix(:,all_indices) \ incident_wave_vector;
      T(out_medium,ci+total_orders) = Tcol(1:length(indices_medium),1);
      T2(out_cube,ci+total_orders) = Tcol((1+length(indices_medium)):length(all_indices),1);

   if verbose
      fprintf(1,'.');
   end
   
   if verbose
       waitbar(ci/ci_max,h)
   end
  
end
end

if verbose
    close(h)
end

if verbose
   fprintf(1,'done!\n');
end

return


