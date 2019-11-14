% convert 
polv=[randn(3,1),1i*randn(3,1)];
thetax=pi*rand(3,1);
lmod=randi(4,3,1)-2;

%%

%% bessel test cases nargin 3
[n,m,a,b]=bsc_bessel(2,thetax,polv);
[nn,mm,aa,bb]=bsc_bessel_old(2,thetax,polv);

a-aa
b-bb

%% bessel test cases nargin 4

[n,m,a,b]=bsc_bessel(2,thetax,polv(:,1),polv(:,2));
[nn,mm,aa,bb]=bsc_bessel_old(2,thetax,polv(:,1),polv(:,2));

a-aa
b-bb

[n,m,a,b]=bsc_bessel(2,thetax,ones(size(thetax)),polv);
[nn,mm,aa,bb]=bsc_bessel_old(2,thetax,polv,ones(size(thetax)));
a-aa
b-bb

%% bessel test cases nargin 5

[n,m,a,b]=bsc_bessel(2,thetax,-1*ones(size(thetax)),polv(:,1),polv(:,2));
[nn,mm,aa,bb]=bsc_bessel_old(2,thetax,polv(:,1),polv(:,2),-1*ones(size(thetax)));

a-aa
b-bb

