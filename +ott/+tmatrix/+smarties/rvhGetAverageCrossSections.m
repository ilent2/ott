function stCrossSection = rvhGetAverageCrossSections(k1, CstTRa)
%% rvhGetAverageCrossSections
% Calculate orientation-averaged cross sections from T-matrix
% 
% pstGetAverageCrossSections(k1,CstTRa) Calculate absorption, scattering and
% extinction cross sections, for the case of orientationally-averaged incidence.
% This is for matrices which make use of rvh symmetry.
%
% Input:
%		 k1:		wavevector [L x 1]
%		 CstTR:     Cell containing structures with T matrices (in rvh-block form).
%                   Should be {L x M} and should contain all m from m=0..N
%
% Output:
%       stCrossSection: A structure with three fields
%           - Cext: [L x 1] wavelength-dependent extinction coefficient
%           - Csca: [L x 1] wavelength-dependent scattering coefficient
%           - Cabs: [L x 1] wavelength-dependent absorption coefficient
%
% Dependency: 
% none

L=size(CstTRa,1);
M=size(CstTRa,2);

Cextsum = zeros(L,1);
Cscasum = zeros(L,1);

for lInd=1:L % Loop on lambda
    if length(CstTRa{lInd,1}.st4MTeo.ind1)+length(CstTRa{lInd,1}.st4MTeo.ind2) +1 ~= M
            warning('SMARTIES:missingm','Warning in rvhGetAverageCrossSections: CstTRa does not seem to contain T-matrices for all m');
    end
    for mInd = 1:size(CstTRa,2)
        m = mInd-1;% we are making an assumption about the layout of CstTR here
        if (CstTRa{1, mInd}.st4MTeo.m~=m || CstTRa{1, mInd}.st4MToe.m~=m)
            warning('SMARTIES:missingm','Warning in rvhGetAverageCrossSections: CstTRa does not seem to contain T-matrices for all m');
        end
        mFactor = 1;
        if m == 0
            mFactor = 0.5;
        end

        % From Eq. 5.107 of Mischenko 2002
        Cextsum(lInd) = Cextsum(lInd) + mFactor*(  ...
            sum(diag(CstTRa{lInd, mInd}.st4MTeo.M11)) ...
            + sum(diag(CstTRa{lInd, mInd}.st4MToe.M11)) ...
            + sum(diag(CstTRa{lInd, mInd}.st4MTeo.M22)) ...
            + sum(diag(CstTRa{lInd, mInd}.st4MToe.M22)));
        % From Eq. 5.141 of Mischenko 2002
        Cscasum(lInd) = Cscasum(lInd) + mFactor*(...
            sum(sum(abs(CstTRa{lInd, mInd}.st4MTeo.M11).^2))+...
            sum(sum(abs(CstTRa{lInd, mInd}.st4MTeo.M12).^2))+...
            sum(sum(abs(CstTRa{lInd, mInd}.st4MTeo.M21).^2))+...
            sum(sum(abs(CstTRa{lInd, mInd}.st4MTeo.M22).^2))+...
            sum(sum(abs(CstTRa{lInd, mInd}.st4MToe.M11).^2))+...
            sum(sum(abs(CstTRa{lInd, mInd}.st4MToe.M12).^2))+...
            sum(sum(abs(CstTRa{lInd, mInd}.st4MToe.M21).^2))+...
            sum(sum(abs(CstTRa{lInd, mInd}.st4MToe.M22).^2)));
    end
end

stCrossSection.Cext = -4*pi./k1^2.*real(Cextsum); % [L x 1]

stCrossSection.Csca = 4*pi./k1^2.*Cscasum; % [L x 1]

stCrossSection.Cabs = stCrossSection.Cext - stCrossSection.Csca; % [L x 1]

end
