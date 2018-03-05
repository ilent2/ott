function [unitarity,symmetry]=check_tmatrix(T,absorbing);
%check_tmatrix.m : Checks to see if the T-matrix is valid. Several
%                   conditions need to be satisfied. Just because this
%                   function indicates you are using a t-matrix doesn't
%                   mean you are using the correct one for your problem!
%                   NYI: absorbing.

if nargin<2
    absorbing = 0;
end


figz=0;
textz=0;
%if ~absorbing
%check power:

S=2*T+eye(size(T));

absS1=sqrt(sum(abs(S).^2,1));
absS2=sqrt(sum(abs(S).^2,2));

if figz
    figure
    plot(absS1,'b');
    hold on
    plot(absS2,'r');
    hold off
    
    figure
    imagesc(abs(S'*S)+1e-15)
end
%end

%check symmetry:
nmax=floor(sqrt(size(T,1)/2))

Tl=T([1:nmax*(nmax+2)],nmax*(nmax+2)+[1:nmax*(nmax+2)]);
Tu=T(nmax*(nmax+2)+[1:nmax*(nmax+2)],[1:nmax*(nmax+2)]);

cip=[1:nmax*(nmax+2)].';
[n,m]=combined_index(cip);
cim=combined_index(n,-m);

[M1,M2]=meshgrid(m,m);

shouldbezeros=Tl(cim,cim)-(-1).^(M1+M2).*Tu(cip,cip).';
if figz
    figure
    plot(sqrt(sum(abs(shouldbezeros).^2,1)));
    hold on
    plot(sqrt(sum(abs(shouldbezeros).^2,2)),'r');
    hold off
end
%report output

if textz
    disp('Analysys')
    disp('--------')
    disp(['Non-absorbing/gain unitarity test: ',num2str(sqrt(sum(sum(abs(S'*S-eye(size(S))).^2)))*100/sqrt(nmax)),'% cumulative power deviation.']);
    disp(['TETM symmetry test: ',num2str(sqrt(sum(sum(abs(shouldbezeros).^2)))/sqrt(sum(sum(abs(Tl).^2)))*100),'% cumulative power deviation.']);
end

unitarity=sqrt(sum(sum(abs(S'*S-eye(size(S))).^2)))/sqrt(nmax);
symmetry=full(sqrt(sum(sum(abs(shouldbezeros).^2))));
