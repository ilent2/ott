function [ modeweights col_modes row_modes ] = paraxial_transformation_matrix( paraxial_order, basis_in, basis_out, normal_mode )
% PARAXIAL_TRANSFORMATION_MATRIX produces paraxial beam mode conversion
% in a particular order.
%
% [modeweights, col_modes, row_modes] = ...
%   PARAXIAL_TRANSFORMATION_MATRIX(degree, basis_in, basis_out) or
% [modeweights, col_modes, row_modes] = ...
%   PARAXIAL_TRANSFORMATION_MATRIX( degree, basis_in, basis_out, normal_mode)
%
% inputs:
%
% degree : paraxial degree of modes e.g. gaussian is 0.
% basis_in : 0 vortex LG, 1 vortex HG, [2,xi] vortex IG.
% basis_out : 0 LG, 1 HG, [2,xi] IG.
% normal_mode : 0 is default vortex->non-vortex. (because of toolbox modes)
%               1 makes the conversion non-vortex->non-vortex.
%
% outputs:
%
% modeweights : weights of the conversion basis_in->basis_out. Format: each
%               row is the corresponding mode of the output basis, each
%               column for the input basis, such that
%               conj. tranpose(basis_1 --> basis_2) = basis_2 --> basis_1.
%               holds for (vortex->non-vortex)' == non-vortex->vortex.
% col_modes :   outputs the LG, HG or IG indices corresponding to each
%               COLUMN of the matrix. [p,l], [m,n], [o,m,p].
% row_modes :   outputs the LG, HG or IG indices corresponding to each ROW
%               of the matrix. [p,l], [m,n], [o,m,p].
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:paraxial_transformation_matrix:move', ...
    'This file will move to ott.utils.paraxial_transformation_matrix');
ott_warning('internal');

if nargin==3
    normal_mode=0;
end

% We need to have complete transformations of LG->HG->IG (non-vortex) we
% need three types of matrix. As it is bothersome to re-arrange the matrix
% to find symmetric and anti-symmetric parts I also compute these as well.
% These computations are recursive to avoid clutter in the toolbox
% directory. The code blocks conatin redundant outputs which used to go to 
% look-up functions.

switch 100*normal_mode+10*basis_in(1)+basis_out(1)
    case 000
        modeweights=genLG2IG(paraxial_order,0);
    case 001
        modeweights=genLG2HG(paraxial_order);
    case 002
        modeweights=genLG2IG(paraxial_order,basis_out(2));
    case 010
        modeweights=genLG2IG(paraxial_order,0)*genLG2vHG(paraxial_order)';
    case 011
        modeweights=genLG2HG(paraxial_order)*genLG2vHG(paraxial_order)';
    case 012
        modeweights=genLG2IG(paraxial_order,basis_out(2))*genLG2vHG(paraxial_order)';
    case 020
        modeweights=genLG2IG(paraxial_order,0)*genLG2vIG(paraxial_order,basis_in(2))';
    case 021
        modeweights=genLG2HG(paraxial_order)*genLG2vIG(paraxial_order,basis_in(2))';
    case 022
        modeweights=genLG2IG(paraxial_order,basis_out(2))*genLG2vIG(paraxial_order,basis_in(2))';
    case 100
        modeweights=eye(paraxial_order+1);
    case 101
        modeweights=genLG2HG(paraxial_order)*genLG2IG(paraxial_order,0)';
    case 102
        modeweights=genLG2IG(paraxial_order,basis_out(2))*genLG2IG(paraxial_order,0)';
    case 110
        modeweights=genLG2IG(paraxial_order,0)*genLG2HG(paraxial_order)';
    case 111
        modeweights=genLG2HG(paraxial_order)*genLG2HG(paraxial_order)';
    case 112
        modeweights=genLG2IG(paraxial_order,basis_out(2))*genLG2HG(paraxial_order)';
    case 120
        modeweights=genLG2IG(paraxial_order,0)*genLG2IG(paraxial_order,basis_in(2))';
    case 121
        modeweights=genLG2HG(paraxial_order)*genLG2IG(paraxial_order,basis_in(2))';
    case 122
        modeweights=genLG2IG(paraxial_order,basis_out(2))*genLG2IG(paraxial_order,basis_in(2))';
    otherwise
        ott_warning('external');
        error('unknown parameter')
end


switch basis_out(1)
    case 0
        i3_out=[];
        i2_out=[-paraxial_order:2:paraxial_order].';
        i1_out=floor((paraxial_order-abs(i2_out))/2);
    case 1
        i3_out=[];
        i2_out=[0:paraxial_order].';
        i1_out=paraxial_order-i2_out;   
    case 2
        i3_out=[0:paraxial_order].';
        i2_out=paraxial_order-2*(floor(i3_out/2));
        i1_out=paraxial_order*ones(paraxial_order+1,1);  
        i3_out=iseven(i3_out);
        
end
row_modes=[i1_out,i2_out,i3_out];


switch basis_in(1)
    case 0
        i3_in=[];
        i2_in=[-paraxial_order:2:paraxial_order].';
        i1_in=floor((paraxial_order-abs(i2_in))/2);
    case 1
        i3_in=[];
        i2_in=[0:paraxial_order].';
        i1_in=paraxial_order-i2_in;   
    case 2
        i3_in=[0:paraxial_order].';
        i2_in=paraxial_order-2*(floor(i3_in/2));
        i1_in=paraxial_order*ones(paraxial_order+1,1);  
        i3_in=iseven(i3_in);

end
col_modes=[i1_in,i2_in,i3_in];

ott_warning('external');

end

function [output,LGlookups,HGlookups]=genLG2HG(order_paraxial);
% genLG2HG.m --- LG->HG conversion matrix.
%
% Usage:
%
% [modewieghts,LGlookups,HGlookups] = genLG2HG(paraxial_order)

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
end


function [output,LGlookups,HGlookups]=genLG2vHG(order_paraxial);
% genLG2vHG.m --- LG->vortex HG conversion matrix.
%
% Usage:
%
% [modewieghts,LGlookups,HGlookups] = genLG2vHG(paraxial_order)

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
outputt=zeros(order_paraxial+1);

%we're going to produce the spinor inverse of what was in beijersbergen
%1993.
outputt(1:floor(order_paraxial/2)+1,:) = [block_to_mirror(:,1:ceil(order_paraxial/2)),fliplr(block_to_mirror.*(-1).^(P.'))]; %
outputt(end-ceil(order_paraxial/2)+1:end,:)=flipud(outputt(1:ceil(order_paraxial/2),:));

outputt(ceil(order_paraxial/2)+1:end,2:2:end)=-outputt(ceil(order_paraxial/2)+1:end,2:2:end);
outputt(2:2:end,:)=-outputt(2:2:end,:);
output(1:floor((order_paraxial+1)/2),:)=1/sqrt(2)*(outputt(1:2:floor((order_paraxial+1)/2)*2,:)+outputt(2:2:floor((order_paraxial+1)/2)*2,:));

if ~rem(order_paraxial,2)
    output(floor((order_paraxial+1)/2)+1,:)=outputt(end,:);
end

output(end-floor((order_paraxial+1)/2)+1:end,:)=fliplr(flipud(output(1:floor((order_paraxial+1)/2),:)));
output=flipud(output);
if nargout>1
    [LGlookups,HGlookups]=meshgrid([0:order_paraxial],[0:order_paraxial]);
end
end

function [modeweights,LGlookups,IGlookups]=genLG2IG(order_paraxial,xi)
% genLG2IG.m --- LG->IG conversion matrix for elipticity xi. A
%		parameter of xi=0 gives the non-vortex LG modes.
%		a parameter of xi=1e100 will give pretty HG modes.
%
% Usage:
%
% [modewieghts,LGlookups,IGlookups] = genLG2IG(paraxial_order,xi)

%first create the upper block... these are the fourier coefficients...
[A_n,B_n]=incecoefficients(order_paraxial,xi);

%prepare the index matrices for this upper block.
p=[floor(order_paraxial/2):-1:0];
l=order_paraxial-2*p;

P=repmat(p,[length(p),1]);
L=repmat(l,[length(l),1]);

Nb=(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));
Na=(1+1*(L==0)).*(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));

NA_n=flipud(Na.*A_n);
NB_n=flipud(Nb.*B_n);

%calculate the other blocks:
if rem(order_paraxial,2)==1
    bigA=[fliplr(NA_n),NA_n];
    bigB=[fliplr(NB_n),-NB_n];
    
    modeweights=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+1]).';
    
else
    bigA=[fliplr(NA_n),NA_n(:,2:end)];
    bigB=[fliplr(NB_n),-NB_n(:,2:end)];
    
    modeweights=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+2]).';
    modeweights(end,:)=[];
end

for ii=1:size(modeweights,1)
    modeweights(ii,:)=modeweights(ii,:)/sqrt(sum(abs(modeweights(ii,:)).^2));
end

%create imaginary matrix:
imat=repmat((-1i).^[0:order_paraxial].',[1,order_paraxial+1]);

modeweights=modeweights.*imat;

if nargout>1
    [LGlookups,IGlookups]=meshgrid([0:order_paraxial],[0:order_paraxial]);
end
end

function [modeweights,LGlookups,IGlookups]=genLG2vIG(order_paraxial,xi)
% genLG2vIG.m --- LG->IG conversion matrix for elipticity xi. A
%		parameter of xi=0 gives the vortex LG modes.
%		a parameter of xi=1e100 will give pretty vortex HG modes.
%
% Usage:
%
% [modewieghts,LGlookups,IGlookups] = genLG2vIG(paraxial_order,xi)

%first create the upper block... these are the fourier coefficients...
[A_n,B_n]=incecoefficients(order_paraxial,xi);

%prepare the index matrices for this upper block.
p=[floor(order_paraxial/2):-1:0];
l=order_paraxial-2*p;

P=repmat(p,[length(p),1]);
L=repmat(l,[length(l),1]);

Nb=(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));
Na=(1+1*(L==0)).*(sqrt(factorial(P+L).*factorial(P)).*(-1).^(P+L+(order_paraxial+L.')/2));

NA_n=flipud(Na.*A_n);
NB_n=flipud(Nb.*B_n);

%calculate the other blocks:
if rem(order_paraxial,2)==1
    bigA=[fliplr(NA_n),NA_n];
    bigB=[fliplr(NB_n),-NB_n];
    
    modeweightst=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+1]).';
    
else
    bigA=[fliplr(NA_n),NA_n(:,2:end)];
    bigB=[fliplr(NB_n),-NB_n(:,2:end)];
    
    modeweightst=reshape([bigA.';bigB.'],[order_paraxial+1,order_paraxial+2]).';
    modeweightst(end,:)=[];
end

for ii=1:size(modeweightst,1)
    
    modeweightst(ii,:)=modeweightst(ii,:)/sqrt(sum(abs(modeweightst(ii,:)).^2));
end

modeweights=zeros(size(modeweightst));

modeweights(1:floor((order_paraxial+1)/2),:)=1/sqrt(2)*(modeweightst(1:2:floor((order_paraxial+1)/2)*2,:)+modeweightst(2:2:floor((order_paraxial+1)/2)*2,:));

if ~rem(order_paraxial,2)
    modeweights(floor((order_paraxial+1)/2)+1,:)=modeweightst(end,:);
end

modeweights(end-floor((order_paraxial+1)/2)+1:end,:)=fliplr(flipud(modeweights(1:floor((order_paraxial+1)/2),:)));

if nargout>1
    [LGlookups,IGlookups]=meshgrid([0:order_paraxial],[0:order_paraxial]);
end
end
