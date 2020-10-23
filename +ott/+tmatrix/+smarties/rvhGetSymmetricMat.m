function CstMaSym = rvhGetSymmetricMat(CstMa, CsMatList)
%% rvhGetSymmetricMat
% Uses symmetry of T-matrix to get lower triangle from upper triangle part
% 
% rvhGetSymmetricMat(CstMa,CsMatList) symmetrises the matrices given in
% CsMatList, and leaves other matrices present in CstMa unchanged.
% All the matrices are returned in CstMaSym. This uses the upper triangular
% matrices to deduce the lower triangular parts from the symmetry relations
% of the T-matrix obtained from Eqs. 5.34 and 5.37 of [Mishchenko 2002],
% i.e: T11 = T11.', T22 = T22.', T12 = -T21.', T21 = -T12.'
% This function assumes rvh symmetry.
%
% Dependency: 
% none

if nargin<2
    CsMatList={'st4MT'};
end
numEntries = numel(CstMa);
CstMaSym = CstMa; % copy the entire variable to non-symmetrized matrices (e.g. R)

for ind = 1:numEntries
    for ii=1:numel(CsMatList)
        
        % indices in the cell array
        [cind1, cind2] = ind2sub(size(CstMa), ind);
        
        % variable name, i.e.: st4MT
        ss=CsMatList{ii};
        
        %% First for eo matrix
        st4M=CstMa{cind1, cind2}.([ss, 'eo']);
        ind1=st4M.ind1;
        ind2=st4M.ind2;
        if ~isempty(ind1) && ~isempty(ind2)
            nOffset12=(ind2(1)-ind1(1)+1)/2;
            nOffset21=1-nOffset12;
            % Keep upper triangular (including diagonal)
            UT=triu(st4M.M11,1); % Upper tri without diagonal
            DD=diag(diag(st4M.M11)); % diagonal only
            st4M.M11=UT+DD+UT.'; % symmetrised matrix

            UT=triu(st4M.M22,1); % Upper tri without diagonal
            DD=diag(diag(st4M.M22)); % diagonal only
            st4M.M22=UT+DD+UT.'; % symmetrised matrix

            UT1=triu(st4M.M12,nOffset12); % Upper tri without diagonal
            UT2=triu(st4M.M21,nOffset21);
            st4M.M12=UT1-UT2.'; % symmetrised matrix

            st4M.M21=-st4M.M12.'; % symmetrised matrix
        end

        CstMaSym{cind1, cind2}.([ss, 'eo'])=st4M;

        %% Next for oe matrix
        st4M=CstMa{cind1, cind2}.([ss, 'oe']);
        ind1=st4M.ind1;
        ind2=st4M.ind2;
        if ~isempty(ind1) && ~isempty(ind2)
            nOffset12=(ind2(1)-ind1(1)+1)/2;
            nOffset21=1-nOffset12;
            % Keep upper triangular (including diagonal)
            UT=triu(st4M.M11,1); % Upper tri without diagonal
            DD=diag(diag(st4M.M11)); % diagonal only
            st4M.M11=UT+DD+UT.'; % symmetrised matrix

            UT=triu(st4M.M22,1); % Upper tri without diagonal
            DD=diag(diag(st4M.M22)); % diagonal only
            st4M.M22=UT+DD+UT.'; % symmetrised matrix

            UT1=triu(st4M.M12,nOffset12); % Upper tri without diagonal
            UT2=triu(st4M.M21,nOffset21);
            st4M.M12=UT1-UT2.'; % symmetrised matrix

            st4M.M21=-st4M.M12.'; % symmetrised matrix
        end

        CstMaSym{cind1, cind2}.([ss, 'oe'])=st4M;


    end
end