function B = change_nmax(A,Nmax)
% change_nmax.m
%
% Resizes a T-matrix or beam coefficients vector to a new Nmax
%
% Usage:
% newT = change_nmax(oldT,Nmax);
%
% A warning is issued if the matrix/vector is being truncated and possibly
% significant values are being discarded.
%
% This file is part of the package Optical tweezers toolbox 1.3
% Copyright 2006-2013 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html


warning_error_level = 1e-6;

size_of_matrix = size(A);

if size_of_matrix(2) == 1
    beam_vector = true;
else
    beam_vector = false;
end

total_orders = combined_index(Nmax,Nmax);

if beam_vector
    
    if total_orders > size_of_matrix(1)
        [row_index,col_index,a] = find(A);
        B = sparse(row_index,col_index,a,total_orders,1);
    else    
        B = A(1:total_orders);
    end
     
else
    
    midpoint = size_of_matrix(1)/2;
    
    A11 = A(1:midpoint,1:midpoint);
    A12 = A(1:midpoint,(midpoint+1):end);
    A21 = A((midpoint+1):end,1:midpoint);
    A22 = A((midpoint+1):end,(midpoint+1):end);
    
    if total_orders > midpoint
    
        [row_index,col_index,a] = find(A11);
        B11 = sparse(row_index,col_index,a,total_orders,total_orders);

        [row_index,col_index,a] = find(A12);
        B12 = sparse(row_index,col_index,a,total_orders,total_orders);

        [row_index,col_index,a] = find(A21);
        B21 = sparse(row_index,col_index,a,total_orders,total_orders);

        [row_index,col_index,a] = find(A22);
        B22 = sparse(row_index,col_index,a,total_orders,total_orders);

    else

        B11 = A11(1:total_orders,1:total_orders);
        B12 = A12(1:total_orders,1:total_orders);
        B21 = A21(1:total_orders,1:total_orders);
        B22 = A22(1:total_orders,1:total_orders);
        
    end
    
    B = [ B11 B12; B21 B22 ];

end
    
magA = full(sum(sum(abs(A).^2)));
magB = full(sum(sum(abs(B).^2)));

apparent_error = abs( magA - magB )/magA;

if apparent_error > warning_error_level
    warning([ 'Apparent error of ' num2str(apparent_error) ]);
end

return
