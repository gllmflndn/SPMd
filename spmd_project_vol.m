function [Av,Cv,Sv] = spmd_project_vol(V,A,dim,hld)
% forms maximium intensity projections from a mapped image
% FORMAT [Av,Cv,Sv] = spmd_project_vol(V,[A,dim,hld])
% V     -  Mapped volume
% A     -  4 x 4 transformation matrix; default is identity.
% dim   -  [i j k] are the three dimensions that defines the output
%          images.   The coordinates in 3-D space of the voxels in this
%          image are assumed to range from 1,1,1 to i,j,k.  Default is
%          the image dimensions.
% hld   -  Hold, interpolation (see spm_sample_vol.m)
%_______________________________________________________________________
%
% Creates maximum intensity projections of a volume.  Akin to
% spm_slice_vol, sampled points are transformed by matrix A to find
% points in voxel space.  
%
% With a default matrix A (identity), this will produce projections
% with voxel resolution; if the voxels are anisotropic, this will
% produced 'squashed' appearance. 
%
% see also spm_mip.m
%
%_______________________________________________________________________
% $Id: spmd_project_vol.m,v 1.1 2005/10/13 03:53:55 nichols Exp $
% Based on spm_project_vol.m 1.3 Tom Nichols 01/03/11

if (nargin<2 | isempty(A)), A = spm_matrix([]); end
if (nargin<3 | isempty(dim)), dim = V.dim(1:3); end
if (nargin<4), hld = 0; end

dim = fix(dim(:)');

% Setup Av Cv Sv
Av = repmat(NaN,dim([1 2]));
Cv = repmat(NaN,dim([1 3]));
Sv = repmat(NaN,dim([2 3]));

% Fill in Ac Cv Sv
for x3=1:dim(3),
  img      = spm_slice_vol(V,A*inv(spm_matrix([0 0 -(x3-1)])),dim([1 2]),hld);
  Av       = max(Av,img);
  Cv(:,x3) = max(Cv(:,x3),max(img,[],2));
  Sv(:,x3) = max(Sv(:,x3),max(img,[],1)');
end;

