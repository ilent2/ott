function B = change_nmax(A,Nmax)
% CHANGE_NMAX Resizes T-matrix or beam vector to a new Nmax
%
% newT = CHANGE_NMAX(oldT,Nmax) creates a new T-matrix, if the new Nmax is
% greater, the new T-matrix is sparse with zeros padding the original.
% If Nmax is smaller, the new T-matrix is a truncated version of the
% original matrix.
%
% newa = CHANGE_NMAX(olda,Nmax) creates a new beam vector with new Nmax.
% The input must not be a combined beam vector, i.e. not [a;b], a and b
% must be resized seperatly.  When the beam vector is expanded, the new
% vector is stored as a sparse matrix.
%
% A warning is issued if the matrix/vector is being truncated and possibly
% significant values are being discarded.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('internal');

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

ott_warning('external');

if apparent_error > warning_error_level
    ott_warning('ott:change_nmax:tol', ...
        ['Apparent error of ' num2str(apparent_error)]);
end

return
