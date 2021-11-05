function X = spmd_sample_vols(V,x,y,z,hold)
% returns voxel values from a collection of memory mapped images
% FORMAT X = spmd_sample_vols(V,x,y,z,hold);
% V      -  is a vector of memory mapped image volumes
% x      -  matrix of x coordinates {pixels}
% y      -  matrix of y coordinates {pixels}
% z      -  matrix of z coordinates {pixels}
% hold   -  sets the interpolation method for the resampling.
%           0          Zero-order hold (nearest neighbour).
%           1          First-order hold (trilinear interpolation).
%           2->127     Higher order Lagrange (polynomial) interpolation using
%                      different holds (second-order upwards).
%          -127 - -1   Different orders of sinc interpolation.
% X      -  output image
%
%_______________________________________________________________________
%
% spmd_sample_vols will return the voxel values from the memory mapped
% volumes indicated by V at coordinates x,y,z.  Values from coordinates
% outside the image are set to zero. x, y and z must be matrices of the
% same dimensions
%
% Based on spm_sample_vol.m  2.1
%
%_______________________________________________________________________
% @(#)spmd_sample_vols.m	1.2 02/05/09


if any(any(diff(cat(1,V.dim),1,1),1)&[1,1,1]) %NB: Bombs for single image
	error('images don''t all have the same dimensions'), end
if any(any(any(diff(cat(3,V.mat),1,3),3)))
        error('images do not all have same orientation & voxel size'), end
 
%- Let spm_sample_vol check validity x,y,z
%-----------------------------------------------------------------------

nV  = length(V);
X   = repmat(0,[size(x) nV]);

for t = 1:nV
  
  X(:,:,t) = spm_sample_vol(V(t),x,y,z,hold);
  
end
