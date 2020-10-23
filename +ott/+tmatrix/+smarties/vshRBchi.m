function chix = vshRBchi(n, x)
  %% vshRBchi
% Calculates the Riccati-Bessel function chi_n(x) = x y_n(x)
% 
%	vshRBchi(n,x) calculates chi_n(x) for many n and x
%	where n is [1 x N] and x is [X x 1]
%   Returns chix [X x N]
%
% Dependency: 
% none

chix = zeros(length(x),length(n));
for xind=1:length(x)
    yx = bessely(n+1/2,x(xind));

    if find(isinf(yx)==true,1) % beyond double precision
        warning(['Warning: Bessel (y) calculation went beyond precision in ', mfilename]);
        warning(['x=', num2str(x(xind)), ' Nmax=', int2str(max(n))]);
    end

    chix(xind,:) = sqrt(x(xind)*pi/2).*yx;
end



end
