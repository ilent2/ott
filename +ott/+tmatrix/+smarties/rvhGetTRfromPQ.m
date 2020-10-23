function CstTRa = rvhGetTRfromPQ(CstPQa, bGetR)
%% rvhGetTRfromPQ
% Calculates T (and possibly R) from P,Q for scatterers with a plane of symmetry
% 
% rvhGetTRfromPQ(CstPQa) Calculates the matrices T,R from the matrices P,Q.
% This makes use of block inversion.
% See Sec. 4.5 and Eq. 70 in JQSRT2013 for details.
% The matrices in CstPQa must make use of rvh symmetry.
% Input:
%   CstPQa: [1 x M] a cell, each element of which contains a structure (one entry
%           for each m), and which makes use of the rvh symmetry
%   bGetR (optional - default=false): boolean. If true R is computed and
%   returned
%
%	Output:
%           CstTRa: A cell of the same size as CstPQa
%                   which contains structures with fields st4MTeo,st4MToe,
%                   st4MReo and st4MRoe
%
% Dependency: 
% isOctave

import ott.tmatrix.smarties.*;

if nargin<2
    bGetR=false;
end

% ignore bad-conditioning warnings
if isOctave()
    currentWarning = warning('off','Octave:nearly-singular-matrix');
else
    currentWarning = warning('off','MATLAB:nearlySingularMatrix');
end

if isstruct(CstPQa)% we have a struct, make it into a cell
    CstPQt{1} = CstPQa;
    CstPQa = CstPQt;
    clear CstPQt
end

% Note that all matrix inversions of a matrix M will use LU decomposition with
% partial pivoting of the columns (not the rows as normally implemented
% in the Matlab LU function).
% This amounts to finding L,U, and P such that L*U = P*M.'
% See built-in function invertLUcol at the end to see how
% the inversion is carried out using LU with column pivoting.


% Loop over m
numEntries = numel(CstPQa);
CstTRa = cell(size(CstPQa));
for ind = 1:numEntries
    % Get blocks of P and Q matrices
    [ind1, ind2] = ind2sub(size(CstPQa), ind);
    stPQ = CstPQa{ind1, ind2};
    % For Peo,Qeo
    P11ee=stPQ.st4MPeo.M11;
    P12eo=stPQ.st4MPeo.M12;
    P21oe=stPQ.st4MPeo.M21;
    P22oo=stPQ.st4MPeo.M22;
    Q11ee=stPQ.st4MQeo.M11;
    Q12eo=stPQ.st4MQeo.M12;
    Q21oe=stPQ.st4MQeo.M21;
    Q22oo=stPQ.st4MQeo.M22;
    % For Poe,Qoe
    P11oo=stPQ.st4MPoe.M11;
    P12oe=stPQ.st4MPoe.M12;
    P21eo=stPQ.st4MPoe.M21;
    P22ee=stPQ.st4MPoe.M22;
    Q11oo=stPQ.st4MQoe.M11;
    Q12oe=stPQ.st4MQoe.M12;
    Q21eo=stPQ.st4MQoe.M21;
    Q22ee=stPQ.st4MQoe.M22;

    Ne=size(Q11ee,1);
    No=size(Q11oo,1);
    ind1eo=stPQ.st4MQeo.ind1;
    ind2eo=stPQ.st4MQeo.ind2;
    ind1oe=stPQ.st4MQoe.ind1;
    ind2oe=stPQ.st4MQoe.ind2;

    m=stPQ.st4MPeo.m;
    
    if (m==0) % m=0 case separately, because 12 and 21 are empty.
        
        st4MRoe.M11 = invertLUcol(Q11oo);
        st4MToe.M11 = - P11oo * st4MRoe.M11;
        st4MRoe.M22 = invertLUcol(Q22ee); 
        st4MToe.M22 = - P22ee * st4MRoe.M22;
        st4MReo.M11 = invertLUcol(Q11ee);  
        st4MTeo.M11 = - P11ee * st4MReo.M11;
        st4MReo.M22 = invertLUcol(Q22oo);   
        st4MTeo.M22 = - P22oo * st4MReo.M22;

        % Off-diagonal blocks
        st4MTeo.M12 = zeros(Ne,No);
        st4MTeo.M21 = zeros(No,Ne);
        st4MToe.M12 = zeros(No,Ne);
        st4MToe.M21 = zeros(Ne,No);

        if bGetR % R-matrix
            st4MReo.M12 = zeros(Ne,No);
            st4MReo.M21 = zeros(No,Ne);
            st4MRoe.M12 = zeros(No,Ne);
            st4MRoe.M21 = zeros(Ne,No);

        end
    else % m <> 0,  do inversion using Eq. 70

        % For Meo matrices
        Q11Inv = invertLUcol(Q11ee); 

        G1=P11ee*Q11Inv;
        G3=P21oe*Q11Inv;
        G5=Q21oe*Q11Inv;
        F2m1=Q22oo-G5*Q12eo;
        F2 = invertLUcol(F2m1); 

        G2=P22oo*F2;
        G4=P12eo*F2;
        G6=Q12eo*F2;

        st4MTeo.M12=G1*G6-G4;
        st4MTeo.M22=G3*G6-G2;
        st4MTeo.M11=-G1-st4MTeo.M12*G5;
        st4MTeo.M21=-G3-st4MTeo.M22*G5;

        if bGetR
            st4MReo.M12=-Q11Inv*G6;
            st4MReo.M22=F2;
            st4MReo.M11=Q11Inv-st4MReo.M12*G5;
            st4MReo.M21=-st4MReo.M22*G5;
        end

        % For Moe matrices
        Q11Inv = invertLUcol(Q11oo); 

        G1=P11oo*Q11Inv;
        G3=P21eo*Q11Inv;
        G5=Q21eo*Q11Inv;
        F2m1=Q22ee-G5*Q12oe;
        F2 = invertLUcol(F2m1); 

        G2=P22ee*F2;
        G4=P12oe*F2;
        G6=Q12oe*F2;

        st4MToe.M12=G1*G6-G4;
        st4MToe.M22=G3*G6-G2;
        st4MToe.M11=-G1-st4MToe.M12*G5;
        st4MToe.M21=-G3-st4MToe.M22*G5;

        if bGetR
            st4MRoe.M12=-Q11Inv*G6;
            st4MRoe.M22=F2;
            st4MRoe.M11=Q11Inv-st4MRoe.M12*G5;
            st4MRoe.M21=-st4MRoe.M22*G5;
        end

    end

    st4MTeo.m=m;
    st4MTeo.ind1=ind1eo;
    st4MTeo.ind2=ind2eo;
    st4MToe.m=m;
    st4MToe.ind1=ind1oe;
    st4MToe.ind2=ind2oe;
    stTR.st4MTeo=st4MTeo;
    stTR.st4MToe=st4MToe;

    if bGetR
        st4MReo.m=m;
        st4MReo.ind1=ind1eo;
        st4MReo.ind2=ind2eo;
        st4MRoe.m=m;
        st4MRoe.ind1=ind1oe;
        st4MRoe.ind2=ind2oe;
        stTR.st4MReo=st4MReo;
        stTR.st4MRoe=st4MRoe;
        stTR.CsMatList={'st4MT', 'st4MR'};
    else
        stTR.CsMatList={'st4MT'};
    end
    CstTRa{ind1, ind2} = stTR;
end

warning(currentWarning); % reactivates bad-conditioning warnings



end


function X = invertLUcol(B)
% This performs the matrix inversion using LU with partial pivoting of columns
% Note that in Matlab, this is equivalent for near-singular matrices to
% X = transpose(transpose(B)\I), but not to X = I/B which does
% pivoting on rows of B
% It is also equivalent to X=inv(B.').', but not to X=inv(B).
% The code below is more explicit and gives the same results in both Matlab and Octave
    
   [L,U,P] = lu(B.'); % LU decomposition with partial pivoting LU = PB.'
    optsLT.LT=true;
    optsUT.UT=true;
    Y=linsolve(L,P,optsLT); % use lower triangular solver
    X=linsolve(U,Y,optsUT).'; % use upper triangular solver

end
