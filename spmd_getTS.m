function [TS,V] = spmd_getTS(P,x,y,z,Sca)
% Return the time series at the position (x,y,z)
% FORMAT [TS,V] = spmd_getTS(P,x,y,z,Sca)
% Input
% P    -  Matrix of filenames -or- Struct array of file handles
% x    -  matrix of x coordinates {pixels}
% y    -  matrix of y coordinates {pixels}
% z    -  matrix of z coordinates {pixels}
% Sca  -  vector of values to multiply images by
% Output
% TS   -  Time series, dimensionality one greater than x,y,z
% V    -  Matrix of filenames
% To speed up, recycle V.
%________________________________________________________________________
% @(#)spmd_getTS.m	1.4 Tom Nichols 02/10/03

%----------------------------- Functions called -------------------------
%    spm_vol
%    spmd_sample_vols
%------------------------------------------------------------------------


% Let spm_sample_vol check validity x,y,z

if isstruct(P)
    V = P;
else
    V = spm_vol(P);
end
nT  = length(V);


if (nargin>4)
    for i=1:nT
        V(i).pinfo(1:2,:) = V(i).pinfo(1:2,:)*Sca(i);
    end
end

TS = spmd_sample_vols(V,x,y,z,0);
TS = reshape(TS,[prod(size(x)) nT])';
TS = squeeze(reshape(TS,[nT size(x)]));
