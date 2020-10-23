function psix = vshRBpsi(n,x)
%% vshRBpsi
% Calculates the Riccati-Bessel function psi_n(x)=x j_n(x)
% 
%	vshRBpsi(n,x) calculates psi_n(x) for many n and x
%	where n is [1 x N] and x is [X x 1]
%   Returns psix [X x N]
%
% Dependency: 
% none

psix = zeros(length(x),length(n));
for xind=1:length(x)
    jx = besselj(n+1/2,x(xind));

    if ~isempty(find(jx==0,1)) % beyond double precision
        warning(['Warning: Bessel (j) calculation went beyond precision in ', mfilename]);
        warning(['x=', num2str(x(xind)), ' Nmax=', int2str(max(n))]);
    end

    psix(xind,:) = sqrt(x(xind)*pi/2).*jx;
end



end
