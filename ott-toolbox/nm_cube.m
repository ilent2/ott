function [ncube,mcube,cicube]=nm_cube(n,m,Nmax)
% nm_cube.m: Calculates which n, m and ci to use for a cube
%
% PACKAGE INFO

%n = 1;
%m = 2;
%Nmax = 20;

m0 = [sort(-[-m:4:Nmax]) (m+4):4:Nmax];

clear nm

l=1;
for k = 1:length(m0)
    
    if isodd(n) & iseven(m0(k))
        n0 = abs(m0(k))+1:2:Nmax;
    elseif isodd(n) & isodd(m0(k))
        n0 = abs(m0(k)):2:Nmax;
    elseif iseven(n) & iseven(m0(k))
        if m0(k) == 0;
            n0 = abs(m0(k))+2:2:Nmax;
        else
        n0 = abs(m0(k)):2:Nmax;
    end
    elseif iseven(n) & isodd(m0(k))
        n0 = abs(m0(k))+1:2:Nmax;
    end
    nm(1,l:l+length(n0)-1) = n0;
    nm(2,l:l+length(n0)-1) = m0(k)*ones(1,length(n0));
    l = l+length(n0);
    
end

nm = nm';
ci = combined_index(nm(:,1),nm(:,2));
ci = sort(ci);
[nn mm] = combined_index(ci); 
all = [nn mm ci];

ncube = all(:,1);
mcube = all(:,2);
cicube = all(:,3);

%max(size(all))/max(all(:,3))