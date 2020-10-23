function [bGetR,Delta,NB,absmvec,bGetSymmetricT,bOutput] = slvGetOptionsFromStruct(stParams, stOptions)
%% slvGetOptionsFromStruct
% Reads optional parameters from struct, else set to default values
% (used in the "solve" functions)
%
% Dependency:
% none

% Calculates R or not?
if isfield(stOptions,'bGetR')
    bGetR = stOptions.bGetR;
else
    bGetR = false; % default is no R
end
% Delta
if isfield(stOptions,'Delta')
    Delta = stOptions.Delta;
else
    Delta = 0; % default is Delta = 0
end
% NB
if isfield(stOptions,'NB')
    NB = stOptions.NB;
else
    NB = 0; % default is NB = 0 (which means that NB will be estimated)
end
% absmvec
if isfield(stOptions,'absmvec')
    absmvec = stOptions.absmvec;
else
    absmvec = (0:1:stParams.N).'; % default is all m
end
% Symmetrization the T-matrix (using the upper part to get the lower)
if isfield(stOptions,'bGetSymmetricT')
    bGetSymmetricT = stOptions.bGetSymmetricT;
else
    bGetSymmetricT = false; % default is false
end
% Display output
if isfield(stOptions,'bOutput')
    bOutput = stOptions.bOutput;
else
    bOutput = true; % default is true (display output is on)
end

end
