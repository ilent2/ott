function [output,LGlookups,HGlookups]=genLG2HG(order_paraxial);
% genLG2HG.m --- LG->HG conversion matrix.
%
% Usage:
%
% [modewieghts,LGlookups,HGlookups] = genLG2HG(paraxial_order)
%
% PACKAGE INFO

n=[0:floor(order_paraxial/2)];
k=n;
[N,K]=meshgrid(n,k);
M=order_paraxial-N;
P=min(N,M);

%Let's start with the troublesome sum. We need to keep convergence... best
%seems to recur with inner index. We only need (in the uppder block) to:
%
%   i=0:        S_i=m!n!/(m-k+i)!/(n-i)!nchoosek(k,i)
%   0<i<n/2:    S_i=-(n-i)/(m-k+i+1)*i/(i+1)S_(i-1)
%
% which seems stable.

%normalisation matrix:
normalisation_matrix=sqrt(factorial(order_paraxial-K).*factorial(K)./factorial(N)./factorial(M)/2^order_paraxial);

summed_matrix=zeros(size(normalisation_matrix));

%special cases
summed_matrix(:,1)=factorial(M(:,1))./factorial(M(:,1)-K(:,1));
summed_matrix(1,:)=ones([1,floor(order_paraxial/2)+1]);

s_i=zeros(length(k)-1,floor(order_paraxial/2)+1);

for jj=1:length(n)-1
    mm=order_paraxial-jj;
    for kk=1:length(k)-1
        s_i(kk,1)=factorial(jj)*factorial(mm)/factorial(mm-kk)/factorial(jj)*nchoosek(kk,0);
        
        for ii=1:kk
            
            s_i(kk,ii+1)=-(jj-ii+1)/(mm-kk+ii)/(ii)*(kk-ii+1)*s_i(kk,ii);
            
        end
        
    end
    summed_matrix(2:end,jj+1)=sum(s_i,2);
end

%we now multiply for upper block!
block_to_mirror=normalisation_matrix.*summed_matrix./factorial(K).*(-1).^(P);

output=zeros(order_paraxial+1);

%we're going to produce the spinor inverse of what was in beijersbergen
%1993.
output(1:floor(order_paraxial/2)+1,:) = [block_to_mirror(:,1:ceil(order_paraxial/2)),fliplr(block_to_mirror.*(-1).^(P.'))]; %
output(end-ceil(order_paraxial/2)+1:end,:)=flipud(output(1:ceil(order_paraxial/2),:));

output(ceil(order_paraxial/2)+1:end,2:2:end)=-output(ceil(order_paraxial/2)+1:end,2:2:end);
output(2:2:end,:)=-output(2:2:end,:);
output=output.*(1i).^repmat([0:order_paraxial].',[1,order_paraxial+1]);

if nargout>1
	[LGlookups,HGlookups]=meshgrid([0:order_paraxial],[0:order_paraxial]);
end
