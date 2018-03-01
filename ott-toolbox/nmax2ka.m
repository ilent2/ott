function ka = nmax2ka(Nmax)
% nmax2ka.m - Finds size parameter ka corresponding to maximum order at
%             which spherical expansions are truncated
%
% Truncation order is given by Nmax = ka + 3 (ka)^(1/3)
%
% PACKAGE INFO

for ii=1:length(Nmax)
    kas = roots([1 (-3*Nmax(ii)) (27+3*Nmax(ii).^2) (-Nmax(ii).^3)]);
    ka(ii) = kas(3);
end

return
