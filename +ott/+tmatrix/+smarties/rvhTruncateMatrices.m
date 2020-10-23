function CstMat = rvhTruncateMatrices(CstMat, nNmax)
%% rvhTruncateMatrices
% Truncate matrices to a maximum number of multipoles N
%
% rvhTruncateMatrices(CstMat, nNmax) Truncate P,Q,R or T matrices in CstMat
% such that the maximum N is now nNmax. This also removes values of m that
% are larger than nNmax.
% This function works on CstMat which makes use of the rvh symmetry.
%
% Input:
%           CstMat: A cell with one entry for each m, which contains
%                   structures with fields M11, M12, M21, M22, ind1, ind2
%                   and m. m must be scalar, ind1 and ind2 should be
%                   vectors containing the valid indices, while the other
%                   entries should be matrices of appropriate sizes.
%                   Fields called CsMatList are ignored (and not altered).
%           nNmax:  The maximum value of N (and m) that is desired in the
%                   output
%
% Output:
%           CstMat: A cell with the same structure as CstMat (input), but
%           where all M11, M12, M21, M22 have been truncated to N=nNmax, and
%           values of m > nNmax are removed
%
% Dependency:
% none

fieldsToTruncate = {'M11', 'M12', 'M21', 'M22'};

keepList = [];
for ii = 1:numel(CstMat)
    fields1 = fieldnames(CstMat{ii});
    for jj = 1:numel(fields1)
        if ~strcmp(fields1{jj}, 'CsMatList') % this can't be truncated
            fields2 = fieldsToTruncate;
            for kk=1:numel(fields2)
                m = CstMat{ii}.(fields1{jj}).('m');
                if m <= nNmax % otherwise we will remove that entry
                    if isempty(find(keepList==ii, 1))
                        keepList = [keepList, ii];
                    end

                    newSize = nNmax - max(1,m)+1;

                    % we are dealing with the matrix with zeros removed
                    [~,ind1Valid] = find(CstMat{ii}.(fields1{jj}).ind1 <= newSize);
                    [~,ind2Valid] = find(CstMat{ii}.(fields1{jj}).ind2 <= newSize);

                    currentMatrix = CstMat{ii}.(fields1{jj}).(fields2{kk});
                    CstMat{ii}.(fields1{jj}).ind1 = CstMat{ii}.(fields1{jj}).ind1(ind1Valid);
                    CstMat{ii}.(fields1{jj}).ind2 = CstMat{ii}.(fields1{jj}).ind2(ind2Valid);
                    switch fields2{kk}
                        case 'M11'
                            CstMat{ii}.(fields1{jj}).(fields2{kk}) = currentMatrix(ind1Valid, ind1Valid);
                        case 'M12'
                            CstMat{ii}.(fields1{jj}).(fields2{kk}) = currentMatrix(ind1Valid, ind2Valid);
                        case 'M21'
                            CstMat{ii}.(fields1{jj}).(fields2{kk}) = currentMatrix(ind2Valid, ind1Valid);
                        case 'M22'
                            CstMat{ii}.(fields1{jj}).(fields2{kk}) = currentMatrix(ind2Valid, ind2Valid);
                    end
                end
            end
        end

    end
end

CstMat = CstMat(keepList);


end
