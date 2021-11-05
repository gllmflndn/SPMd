function PstRes = spmd_MkResid(xSPM,Dir)
% Compute the studentized residuals
% FORMAT PsRes = spmd_MkResid(xSPM, Dir)
% Input
%  Dir    - Directory of SPM results
%  xSPM   - Parameters from SPM.mat
% Output
%  PstRes - Matrix of filename of (internally) studentized residuals
% _________________________________________________________________________
%
% The output of this function is the matrix of filename of internally
% studentized residuals. On the other hand, it creates the studentized
% residuals images in the data directory which will be use in the model
% detail (spatial detail) window.
%
% In the calculation of the studentized residuals, drift terms (DCTs) are
% also accounts.
% _________________________________________________________________________
%
% This function is modified from MkResid.m by Tom Nichols. The only
% difference is it only creates internally studentized residuals.
%
% _________________________________________________________________________
% @(#)spmd_MkResid.m	1.1 04/07/08

% ______________________ Functions called _________________________________
%
%           spm_select
%           spm_vol
%           spm_type
%           spm_str_manip
%           spm_create_image
%           spm_resss
%           spm_slice_vol
%           spm_write_plane
% _________________________________________________________________________

if (nargin<1)
    swd  = spm_str_manip(spm_select(1,'SPM.mat','Select SPM.mat'),'H');
    xSPM = load(fullfile(swd,'SPM.mat'));
    xSPM = xSPM.SPM;
end

if (nargin<2)
    Dir = spm_select(1,'dir','Select directory to save the residual images!');
end

if strcmp(spm_str_manip(Dir,'t'),'.')
    Dir = spm_str_manip(Dir,'h');
else
    Dir = Dir;
end

fprintf('Configuring...')

xGX    = xSPM.xGX;
xX     = xSPM.xX;
VY     = xSPM.xY.VY;
Vbeta  = xSPM.Vbeta;
VResMS = xSPM.VResMS;

nScan  = length(VY);
nBeta  = length(Vbeta);

%- Adding full name to Vbeta file
for i = 1:size(Vbeta,2)
    if isempty(fileparts(deblank(Vbeta(i).fname)))
        Vbeta(i).fname = fullfile(xSPM.swd,Vbeta(i).fname);
    end
end

%- Adding full name to VResMS file
if isempty(fileparts(deblank(VResMS.fname)))
    VResMS.fname = fullfile(xSPM.swd,VResMS.fname);
end

VVar   = spm_vol(VResMS);

X = xX.X;
W = xX.W;


%- Account for filter
%-------------------------------------------------------
DesMtx = xX.X;
if isstruct(xX.K)
    K = xX.K;
    nPS = 0;     % Number of scans in previous sessions
    
    for s=1:length(xX.K)
        KH = xX.K(s).X0;
        nSS = size(KH,1); % Number of scans in this session
        nH  = size(KH,2); % Number of HP bases
        
        %- Zero pad to account for other sessions
        %--------------------------------------------------
        KH = [zeros(nPS,nH); full(KH); zeros(nScan-nPS-nSS,nH)];
        DesMtx = [DesMtx KH];
        nPS = nPS + nSS;
    end
    
else
    K = [];
end

VRs = VY;

for i=1:nScan
    
    %- Impliment whatever scaling was done in spm_spm.m
    %---------------------------------------------------
    VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*xGX.gSF(i);
    
    
    %- Initialize standardized residual images
    %---------------------------------------------------
    [pth,nm,xt,vr] = fileparts(deblank(VY(i).fname));
    VRs(i).fname = fullfile(Dir,['e' nm xt vr]);
    %  We trust that 'VRs(i).n' is set appropriately
    VRs(i).dt = [spm_type('float32') spm_platform('bigend')];
    VRs(i).pinfo = [1 0 0]';
    VRs(i) = spm_create_vol(VRs(i));
    
end

%- Make raw residuals (for now, write into files of stdzd res)
%-----------------------------------------------------
H = W*DesMtx*pinv(W*DesMtx);
I = eye(length(W));
VRs = spmd_resss(VY,VRs,(I-H)*W,'');

%- Make (internally) standardized residuals
%-----------------------------------------------------
fprintf('Creating (internally)  standardized residual... ')
hii = diag(H);
for z=1:VVar.dim(3)
    fprintf('%s%16s',repmat(sprintf('\b'),1,16),...
        sprintf('...plane %3d/%-3d',z,VVar.dim(3)))       %-#
    Var       = spm_slice_vol(VVar,spm_matrix([0 0 z]),VVar.dim(1:2),0);
    for i=1:nScan
        % Read residuals
        Res     = spm_slice_vol(VRs(i),spm_matrix([0 0 z]),VRs(i).dim(1:2),0);
        % Only standardize if residuals has non-zero variance (i.e. exlude
        % 'exact-zero residuals arising from perfect fit).
        if (1-hii(i))>sqrt(eps)
            % Studentized residual
            SRes    = Res./sqrt(Var*(1-hii(i)));
        else
            % Res should be zero, but just in case.
            SRes = Res./sqrt(Var);
        end
        % Write out
        VRs(i)  = spm_write_plane(VRs(i),SRes,z);
    end
end;
%VRs         = spm_close_vol(VRs);
for i = 1:length(VRs)
    PstRes{i} = [VRs(i).fname ',' num2str(VRs(i).n(1))];
end
PstRes      = strvcat(PstRes);
