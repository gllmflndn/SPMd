function spmd_check_registration(images,tag)
% A visual check of image registration quality.
% FORMAT spmd_check_registration or spmd_check_registration(images)
%        Orthogonal views of one or more images are displayed.  Clicking in
%        any image moves the centre of the orthogonal views.  Images are
%        shown in orientations relative to that of the first selected image.
%        
%        The first specified image is shown at the top-left, and the last at
%        the bottom right.  The fastest increment is in the left-to-right
%        direction (the same as you are reading this).
% FORMAT spmd_check_registration(images,Tag)
%  images  - matrix of the names of the images.
%  Tag     - Tag of the window for the displayed images.
%_______________________________________________________________________
% Notices:
%  This program is modified from the spm_check_registration.m by John
%  Ashburner for the diagnostic used. Using this function, multiple
%  windows of the images could be created and tagged "Tag". 
%_______________________________________________________________________
% @(#)spm_check_registration.m  2.2 John Ashburner 99/10/29
% @(#)spmd_check_registration.m 1.8 Wen-Lin Luo 02/10/03


if ~nargin
    images = spm_select([1 15],'image','Select images');
    spmd_check_registration(images);
    return;
end

if nargin < 2, tag = 'Graphics'; end

fg = spm_figure('GetWin',tag);
spm_figure('Clear',fg);
spmd_orthviews('fig',fg,'Reset');
if nargin > 1, spmd_SptlBar; end

mn = size(images,1);
n  = round(mn^0.4);
m  = ceil(mn/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;
for ij=1:mn
    i  = 1-h*(floor((ij-1)/n)+1);
    j  = w*rem(ij-1,n);
    spmd_orthviews('fig',fg,'Image',deblank(images(ij,:)),[j+ds/2 i+ds/2 w-ds h-ds]);
    if ij==1, spmd_orthviews('fig',fg,'Space'); end
end
