function D=half_integer_wigner_rotation_matrix(Nmax,R);
%Generates the orthogonal half-integer order Wigner D matrix for rotation
%matrix R. It does this by a recursion found in Quantum Theory of Angular
%Momentum by Varshalovich et al. pp. 92.

D=spalloc(Nmax*(Nmax+1),Nmax*(Nmax+1),2*Nmax);
t=acos((trace(R)-1)/2);

n=[(R(3,2)-R(2,3))/2,(R(1,3)-R(3,1))/2,(R(2,1)-R(1,2))/2];
if t==0
    n=[0,0,0];
end

[alpha,beta,gamma]=axisangletoeuler(n,t);

iapgon2=1i*(alpha+gamma)/2;
iamgon2=1i*(alpha-gamma)/2;

%the first element is P_0(cos\beta) = 1 so the first block in D is
%for m'!=n
D(1:2,1)=[cos(beta/2)*exp(-iapgon2);-sin(beta/2)*exp(-iamgon2)];
%for m'!=-n
D(1:2,2)=[sin(beta/2)*exp(iamgon2);cos(beta/2)*exp(iapgon2)];
%which is the 2d spinor rotation!.

W=wigner_rotation_matrix(Nmax,R);

for ii=2:Nmax
    
    elementsforW=(ii-1)^2+[0:2*(ii-1)];
    elementsforD=ii*(ii-1)+[1:2*ii];
    
    %There are four different expressions. The J=M, J=-M, J=M', J=-M' are
    %all special cases. We therefore need four different expressions.
    
    %First create a tempory container for all the M, M' we have:
    Dtemp=zeros(2*ii);
    
    %We'll use the two recusrions found in Quantum Theory of Angular
    %Momentum pp. 92, (14) and (15).
    
    %J-M or J+M is simply the count down or up to 2*J from zero.
    JpM=repmat(sqrt([1:2*ii-1]),[2*ii-1,1]);
    
    % We only need a different recursion for either J=-M' or J=M' We'll
    % choose M'... Note I'ver made a transcription error sometwhere and
    % have the indexes muddled up along the diagonalz. I swapped the
    % exp(iapgon2) and exp(-iapgon2) terms.
    Dtemp(1:2*ii-1,1:2*ii-1)=flipud(JpM')./fliplr(JpM).*W(elementsforW,elementsforW)*cos(beta/2)*exp(-iapgon2);
    Dtemp(2:2*ii,1:2*ii-1)=Dtemp(2:2*ii,1:2*ii-1)-JpM'./fliplr(JpM).*W(elementsforW,elementsforW)*sin(beta/2)*exp(-iamgon2);
    
    %last col
    Dtemp(1:2*ii-1,2*ii)=flipud(JpM(2*ii-1,1:2*ii-1)')./JpM(1:2*ii-1,2*ii-1).*W(elementsforW,elementsforW(end))*sin(beta/2)*exp(iamgon2);
    Dtemp(2:2*ii,2*ii)=Dtemp(2:2*ii,2*ii)+JpM(2*ii-1,1:2*ii-1)'./JpM(1:2*ii-1,2*ii-1).*W(elementsforW,elementsforW(end))*cos(beta/2)*exp(iapgon2);
    
    D(elementsforD,elementsforD)=Dtemp;
end

return