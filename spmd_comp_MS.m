function spmd_comp_MS(SPM,varargin)
% Compute diagnosis statistics for general linear model
% FORMAT spmd_comp_MS
% FORMAT spmd_comp_MS(xSPM,varargin)
%  xSPM     - all parameters from SPM.mat
%  varargin - cell arrary of statistics' name
%__________________________________________________________________________
%
% There are three assumptions for general linear model: independence,
% constant variance, and normality. The diagnosis statistics computed here
% test these three assumptions at every voxel and thus create the diagnosis
% images for spatial summary viewer in the SPMd toolbox.
%
% spmd_comp_MS is the main function to compute all diagnosis statistics for
% general linear model in the spatial summary viewer. Before using this
% function, output files, SPM.mat, from SPM should exist in the current
% directory for the calculation of statistics.
%__________________________________________________________________________
%
% The output of spmd_comp_MS takes the form of an SPMd_(Diagnosis).mat of
% diagnosis statistics, and diagnostic images of the statistics and
% corresponding log10(P-value).
%
% SPMd_(Diagnosis).mat files SPMd_P(Diagnosis).mat - diagnosis statistics
% SPMd_Corr SPMd_Pcorr
%            - cellstr of image filename for Durbin-Waston statistic
% SPMd_Dep  SPMd_PDep
%            - cellstr of image filename for Cumulative periodogram
%              statistic
% SPMd_Homo1 SPMd_PHomo1
%            - cellstr of image filename for Cook-Weisberg statistic
%              score statistics: With respect to experimental design
% SPMd_Homo2 SPMd_PHomo2
%            - cellstr of image filename for Cook_Weisberg statistic
%              score statistics: With respect to predicted response
% SPMd_Homo3 SPMd_PHomo3
%            - cellstr of image filename for Cook_Weisberg statistic
%              score statistics: With respect to global signal
% SPMd_Norm  SPMd_PNorm
%            - cellstr of image filename for Shapiro-Wilk statistic
% SPMd_Outl  SPMd_Poutl
%            - cellstr of image filename for Spatial outlier proportion
%
%__________________________________________________________________________
%
%
% SPMd_Corr.nii SPMd_PCorr.nii    -Durbin_Watson statistic
% These are images of the Durbin-Watson statistic and corresponding
% log10(P-value). This statistic is used to test the assumption of
% independence. Voxels outside the analysis mask are given value NaN.
%
% SPMd_Dep.nii SPMd_PDep.nii      -Cumulative periodogram statistic
% These are images of the cumulative periodogram statistic and
% corresponding log10(P-value). This statistic is used to test the
% non-AR(1) correlation structure. Voxels outside the analysis Mask are
% given NaN.
%
% SPMd_Homo1.nii SPMd_PHomo1.nii  -Cook-Weisberg score statistic
% These are images of the Cook-Weisberg score statistic with respect to
% experimental matrix and corresponding log10(P-value). This is used to
% test if the residuals depend on design matrix. Voxels outside the
% analysis mask are given value NaN.
%
% SPMd_Homo2.nii SPMd_PHomo2.nii  -Cook-Weisberg score statistic
% These are images of the Cook-Weisberg score statistic with respect to
% predicted response and corresponding log10(P-value). This is used to test
% if the residuals depend on predicted response. Voxels outside the
% analysis mask are given value NaN.
%
% SPMd_Homo3.nii SPMd_PHomo3.nii  -Cook-Weisberg score statistic
% These are 3images of the Cook-Weisberg score statistic with respect to
% global signal and corresponding log10(P-value). This is used to test if
% the residuals depend on global signal. Voxels outside the analysis mask
% are given value NaN.
%
% SPMd_Norm.nii SPMd_PNorm.nii    -Shapiro-Wilk statistic
% These are images of the Shapiro-Wilk statistic and corresponding
% log10(P-value). This statistic is used to test the assumption of
% normality. Voxels outside the analysis mask are given value NaN.
%
% SPMd_Outl.nii SPMd_POutl.nii    -Spatial outlier proportion
% These are images of the spatial outlier proportion and corresponding
% log10(P-value). Voxels outside the analysis mask are given value NaN.
%
% SPMd_ResRMS.nii SPMd_PResRMS.nii -sqrt(MSE)
% These are images of the spatial square root of sum of squared mean
% residuals and corresponding log10(P-value). Voxels outside the analysis
% mask are given value NaN.
%__________________________________________________________________________
%
% Reference:
% Luo, W-L and Nichols T. E. (2002) Diagnosis and Exploration of
% Massively Univariate fMRI Models. NeuroImage,19:1014-1032, 2003
%__________________________________________________________________________


%-Say hello
%-----------------------------------------------------------------------
Finter   = spm('FigName','Stats: estimation...');

%-Get SPM.mat if necessary
%-----------------------------------------------------------------------
if nargin == 0
    swd  = spm_str_manip(spm_select(1,'SPM.mat','Select SPM.mat'),'H');
    load(fullfile(swd,'SPM.mat'));
else
    swd     = pwd;
end

%-Ensure data are assigned
%-----------------------------------------------------------------------
try
    SPM.xY.VY;
catch
    helpdlg({   'Please assign data to this design',...
        'Use fMRI under model specification'});
    spm('FigName','Stats: done',Finter);
    return
end

%- check whether the data is multi-session or not?
%- for group data analysis, there is no session parameter,
%- we have to check this property
SPMnames = fieldnames(SPM);
flagSPMnames = 0;
for i=1:length(SPMnames)
    if strcmp(SPMnames(i),'Sess')
        flagSPMnames = 1;
    end
end

if flagSPMnames
    if length(SPM.Sess) > 1
        helpdlg({'OOPs! SPMd can not analysis multi-session data currently'});
        spm('FigName','Stats: done',Finter);
        return
    end
end

%-Delete files from previous analyses
%-----------------------------------------------------------------------
if exist(fullfile('.','SPMd.mat'),'file') == 2
    str   = {'Current directory contains SPM diagnostic files:',...
        'pwd = ',pwd,...
        'Existing results will be appended!'};
    abort = spm_input(str,1,'bd','stop|continue',[1,0],1);
    if abort
        spm('FigName','Stats: done',Finter);
        return
    else
        str = sprintf('Appending old results\n\t (pwd = %s) ',pwd);
        warning(str)
        drawnow
    end
end

%-Specify all the spatial diagnostic statistics
%-----------------------------------------------------------------------
D = get_diag_stat(SPM);

%-Specify the action for calculation of diagnostic statistics
%-----------------------------------------------------------------------
Id   = [];
s    = 0;
if length(varargin) == 0
    stat   = D;
    aPower = 0;
    Tout   = zeros(SPM.nscan,1);
else
    for i = 1:length(varargin)
        Im = varargin{i};
        for j = 1:length(Im)
            Id(s+j)    = strmatch(Im(j),{D(:).name});
            if strcmp(D(Id(s+j)).name,'Dep')
                aPower = 0;
            elseif strcmp(D(Id(s+j)).name,'Outl')
                Tout   = zeros(SPM.nscan,1);
            end
        end
        s  = s + length(Im);
    end
    stat = D(Id);
end
nStat = size(stat,2);

%=======================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%=======================================================================

%-Initialise
%=======================================================================
fprintf('%-40s: %30s','Initialising parameters','...computing')    %-#
xX             = SPM.xX;
df             = SPM.xX.erdf;
[nScan, nBeta] = size(xX.X);

%-If xM is not a structure then assumme it's a vector of thresholds
%-----------------------------------------------------------------------
try
    xM = SPM.xM;
catch
    xM = -ones(nScan,1)/0;
end
if ~isstruct(xM)
    xM = struct('T',    [],...
        'TH',   xM,...
        'I',    0,...
        'VM',   {[]},...
        'xs',   struct('Masking','analysis threshold'));
end

%-Check all necessary variables
%-----------------------------------------------------------------------
xVi   = SPM.xVi;  %-covariance components
V     = xVi.V;    %-Get non-sphericity V
W     = xX.W;     %-Get whitening/Weighting matrix

%-Design space and projector matrix [pseudoinverse] for WLS
%=======================================================================
xX.xKXs = spm_sp('Set',spm_filter(xX.K,W*xX.X));    % KWX
xX.pKX  = spm_sp('x-',xX.xKXs);             % projector


%-Image dimensions and data
%=======================================================================
VY       = SPM.xY.VY;
M        = VY(1).mat;
DIM      = VY(1).dim(1:3);
xdim     = DIM(1); ydim = DIM(2); zdim = DIM(3);
YNaNrep  = VY(1).dt(2);

%-Maximum number of residual images for smoothness estimation
%-----------------------------------------------------------------------
MAXRES   = spm_get_defaults('stats.maxres');

%-maxMem is the maximum amount of data processed at a time (bytes)
%-----------------------------------------------------------------------
MAXMEM   = spm_get_defaults('stats.maxmem');
nSres    = min(nScan,MAXRES);
blksz    = ceil(MAXMEM/8/nScan);                %-block size
nbch     = ceil(xdim*ydim/blksz);               %-# blocks

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')             %-#

%-Initialise output images (unless this is a 1st pass for ReML)
%=======================================================================
fprintf('%-40s: %30s','Output images','...initialising')             %-#

%-Intialise the name of the mean image
%-----------------------------------------------------------------------
VMean = struct(...
    'fname',   'Mean.nii',...
    'dim',     DIM,...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spmd: Mean');

VMean = spm_data_hdr_write(VMean);

%-Intialise standard deviation of residual image file
%-----------------------------------------------------------------------
VResRMS = struct(...
    'fname',   'SPMd_ResRMS.nii',...
    'dim',     DIM,...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spmd: ResRMS');

VResRMS = spm_data_hdr_write(VResRMS);

%-Intialise diagnostic statistics image files
%-----------------------------------------------------------------------
Vstat(1:nStat) = deal(...
    struct(...
    'fname',   '',...
    'dim',     DIM,...
    'dt',      [spm_type('float64'), spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', ''));

for i = 1:nStat
    Vstat(i).fname   = sprintf('SPMd_%s.nii',stat(i).name);
    Vstat(i).descrip = sprintf('SPM_Diag: %s',stat(i).help);
    spm_unlink(Vstat(i).fname)
end

Vstat = spm_data_hdr_write(Vstat);

%-Initialise diagnostic statistics corresponding P-value image files
%-----------------------------------------------------------------------
%VPstat(1:nStat) = deal(...
%    struct(...
%    'fname',   [],...
%    'dim', [DIM',spm_type('float')],...
%    'mat', M,...
%    'pinfo',   [1 0 0]',...
%    'descrip', ''));

VPstat(1:nStat) = deal(...
    struct(...
    'fname',   '',...
    'dim',     DIM,...
    'dt',      [spm_type('float64'), spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', ''));

for i = 1:nStat
    VPstat(i).fname   = sprintf('SPMd_P%s.nii',stat(i).name);
    VPstat(i).descrip = sprintf('SPM_Diag: %s',stat(i).Pdesp);
    spm_unlink(VPstat(i).fname)
end

VPstat = spm_data_hdr_write(VPstat);

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised')      %-#

%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================

%-Intialise variables used in the loop
%=======================================================================
xords = [1:xdim]'*ones(1,ydim); xords = xords(:)';  % plane X coordinates
yords = ones(xdim,1)*[1:ydim];  yords = yords(:)';  % plane Y coordinates
S     = 0;                                          % Volume (voxels)

%-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
%-----------------------------------------------------------------------
XYZ   = zeros(3,xdim*ydim*zdim);

%-Cycle over bunches blocks within planes to avoid memory problems
%=======================================================================
str = 'Diagnostic statistics computation';
spm_progress_bar('Init',100,str,'');

for z = 1:zdim              %-loop over planes (2D or 3D data)
    
    % current plane-specific parameters
    %-------------------------------------------------------------------
    zords   = z*ones(xdim*ydim,1)'; %-plane Z coordinates
    CrBl    = [];                       %-parameter estimates
    CrResI  = [];                       %-normalized residuals
    CrResSS = [];                       %-residual sum of Squares
    
    for i = 1:nStat
        stat(i).Cr = [];                %-diagnostic statistic(s)
        stat(i).P  = [];                %-diagnostic statistic(s)
    end
    Q            = [];          %-in mask indices for this plane
    CrMean       = [];          %-parameter estimates
    
    for bch = 1:nbch            %-loop over blocks
        
        %-# Print progress information in command window
        %---------------------------------------------------------------
        str   = sprintf('Plane %3d/%-3d, block %3d/%-3d',z,zdim,bch,nbch);
        fprintf('\r %-40s: %30s',str,' ')
        
        %-construct list of voxels in this block
        %---------------------------------------------------------------
        I     = [1:blksz] + (bch - 1)*blksz;            %-voxel indices
        I     = I(I <= xdim*ydim);          %-truncate
        xyz   = [xords(I); yords(I); zords(I)];         %-voxel coordinates
        nVox  = size(xyz,2);                    %-number of voxels
        
        %-Get data & construct analysis mask
        %===============================================================
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...read & mask data ')
        Cm    = true(1,nVox);           %-current mask
        
        %-Compute explicit mask
        % (note that these may not have same orientations)
        %---------------------------------------------------------------
        for i = 1:length(xM.VM)
            
            %-Coordinates in mask image
            %-------------------------------------------------------
            j      = xM.VM(i).mat\M*[xyz;ones(1,nVox)];
            
            %-Load mask image within current mask & update mask
            %-------------------------------------------------------
            Cm(Cm) = spm_get_data(xM.VM(i),j(:,Cm)) > 0;
        end
        
        %-Get the data in mask, compute threshold & implicit masks
        %---------------------------------------------------------------
        Y     = zeros(nScan,nVox);
        for i = 1:nScan
            
            %-Load data in mask
            %-------------------------------------------------------
            if ~any(Cm), break, end         %-Break if empty mask
            Y(i,Cm)  = spm_get_data(VY(i),xyz(:,Cm));
            
            Cm(Cm)   = Y(i,Cm) > xM.TH(i);      %-Threshold (& NaN) mask
            if xM.I && ~YNaNrep && xM.TH(i) < 0         %-Use implicit mask
                Cm(Cm) = abs(Y(i,Cm)) > eps;
            end
        end
        
        %-Mask out voxels where data is constant
        %---------------------------------------------------------------
        Cm(Cm) = any(diff(Y(:,Cm),1));
        Y      = Y(:,Cm);               %-Data within mask
        CrS    = sum(Cm);               %-# current voxels
        
        %===============================================================
        %-Proceed with General Linear Model (if there are voxels)
        %===============================================================
        if CrS
            
            %-Whiten/Weight data and remove filter confounds
            %-------------------------------------------------------
            fprintf('%s%30s',repmat(sprintf('\b'),1,30),'filtering')  %-#
            
            KWY   = spm_filter(xX.K,W*Y);
            
            %-General linear model: Weighted least squares estimation
            %------------------------------------------------------
            fprintf('%s%30s',repmat(sprintf('\b'),1,30),'estimation') %-#
            
            beta   = xX.pKX*KWY;        %-Parameter estimates
            res    = spm_sp('r',xX.xKXs,KWY);   %-Residuals
            ResSS  = sum(res.^2);               %-Residual SSQ
            PredY  = KWY - res;             %-Predicted responses
            Mean   = mean(Y);               %-Mean responses
            clear KWY                           %-Clear to save memory
            
            
            %-if we are saving the WLS parmameters
            %------------------------------------------------------
            if isfield(xX,'W')
                
                %-sample covariance and mean of Y (all voxels)
                %----------------------------------------------------
                CrResSS = [CrResSS, ResSS];
            end % (xX,'W')
            clear Y;                            %-Clear to save memory
            
            %-General linear model: Diagnostic processes
            %------------------------------------------------------
            for i = 1:nStat
                if strcmp(stat(i).name,'Homo2')
                    stat(i).X = PredY;
                    [StatD,p,coef] = spmd_stat(stat(i).diag,res,stat(i).X, ...
                        stat(i).para,df);
                elseif strcmp(stat(i).name,'Outl')
                    [StatD,p,coef] = spmd_stat(stat(i).diag,res,stat(i).X, ...
                        stat(i).para,df);
                    Tout           = Tout + StatD{2};
                elseif strcmp(stat(i).name,'Dep')
                    [StatD,p,coef] = spmd_stat(stat(i).diag,res,stat(i).X, ...
                        stat(i).para,df);
                    aPower         = aPower + StatD{2};
                    freq           = StatD{3};
                else
                    [StatD,p,coef] = spmd_stat(stat(i).diag,res,stat(i).X, ...
                        stat(i).para,df);
                end
                
                Diag(i,:)      = StatD{1};
                stat(i).para   = coef;
                Pvalue(i,:)    = p;
            end
            
            %-Save diagnostic statistics for current plane as we go along
            %-----------------------------------------------
            for i = 1:nStat
                stat(i).Cr = [stat(i).Cr, Diag(i,:)];
                stat(i).P  = [stat(i).P, Pvalue(i,:)];
            end
            CrMean       = [CrMean, Mean];
            clear Diag Pvalue
            
        end % (CrS)
        
        %-Append new inmask voxel locations and volumes
        %---------------------------------------------------------------
        XYZ(:,S + [1:CrS]) = xyz(:,Cm);     %-InMask XYZ voxel coords
        Q                  = [Q I(Cm)];     %-InMask XYZ voxel indices
        S                  = S + CrS;       %-Volume analysed (voxels)
        
    end   % (bch)
    
    %-Plane complete, write plane to image files (unless 1st pass)
    %===================================================================
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...saving plane')    %-#
    
    j   = NaN*ones(xdim,ydim);
    
    %-Write diagnostic mean statistics images
    %-------------------------------------------------------------------
    if length(Q), j(Q) = CrMean; end
    VMean   = spm_write_plane(VMean, j, z);
    
    %-Write diagnostic ResRMS statistics images
    if length(Q), j(Q) = sqrt(CrResSS); end
    VResRMS = spm_write_plane(VResRMS,j, z);
    
    
    %-Write diagnostic statistics images
    %-------------------------------------------------------------------
    for i = 1:nStat
        if length(Q), j(Q) = stat(i).Cr; end
        Vstat(i) = spm_write_plane(Vstat(i), j, z);
        if length(Q), j(Q) = stat(i).P; end
        VPstat(i) = spm_write_plane(VPstat(i), j, z);
    end
    
    %-Report progress
    %-------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done')           %-#
    spm_progress_bar('Set',100*(bch + nbch*(z - 1))/(nbch*zdim));
    
end % (for z = 1:zdim)

fprintf('\n')                                                        %-#

for i = 1:nStat
    stat(i).Cr  = [];
    stat(i).P   = [];
    if strcmp(stat(i).name,'Homo2')
        stat(i).X = [];
    end
end
spm_progress_bar('Clear')

%=======================================================================
% - P O S T   E S T I M A T I O N   C L E A N U P
%=======================================================================
if S == 0, warning('No inmask voxels - empty analysis!'), end

if ~exist('Tout')
    Tout = [];
end

if ~exist('aPower')
    aPower = [];
    freq   = [];
end

if isempty(D(5).X)
    Glob   = [];
else
    Glob   = D(5).X(:,2);
end

%-close written image files (unless 1st pass)
%=======================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...closing files')    %-#

%-Set VResRMS scalefactor as sqrt(1/trRV) (raw voxel data is sqrt(ResSS))
%-----------------------------------------------------------------------
VResRMS.pinfo(1) = sqrt(1/xX.trRV);
VResRMS          = spm_data_hdr_write(VResRMS);

%-Save remaining results files and analysis parameters
%=======================================================================
fprintf('%-40s: %30s','Saving results','...writing')                 %-#

%-place fields in SPMd
%-----------------------------------------------------------------------
%-Model summary

if exist(fullfile('.','SPMd_MS.mat'),'file') == 2
    load SPMd_MS
    Istat   = [];
    for j = 1:length(MS.stat)
        Istat(j) = sum(strcmp(MS.stat(j).name,{stat(:).name}));
    end
    MS.stat(Istat==1)   = [];
    MS.Vstat(Istat==1)  = [];
    MS.VPstat(Istat==1) = [];
    nstat = length(MS.stat);
else
    nstat = 0;
    xX.S     = S;           %-Volume (voxels)
    xX.fDM   = D(7).X;              %-Full model design Mtx including
    %intercept and high pass filter
    MS.xX    = xX;              %-design structure
    MS.xVi   = xVi;             %-non-sphericity structure
    MS.VY    = VY;
    MS.VMean    = VMean;        %-Filehandle - Mean
    
    MS.SPM.desp  = 'SPM.mat file used for diagnostic analysis';
    MS.SPM.swd   = swd;
end

for i = 1:nStat
    MS.stat(nstat+i)   = stat(i);
    MS.Vstat(nstat+i)  = Vstat(i);      %-Filehandle - statistics
    MS.VPstat(nstat+i) = VPstat(i);         %-Filehandle - P-value
end

%-Scan Summary
if exist(fullfile('.','SPMd_MS.mat'),'file') == 2
    load SPMd_MS
end
PG.power     = aPower/S;            %-Temporal mean periodogram
PG.freq      = freq;                    %-Frequency of periodogram
GX.ts        = Glob;                %-Global signal
Toutl.ts     = Tout;
xX.fDM       = D(7).X;
SS.xX        = xX;
SS.VY        = VY;
SS.PG        = PG;
SS.GX        = GX;
SS.Toutl     = Toutl;               %-Temporal outlier counts
SS.SPM.swd   = swd;

%-Save analysis parameters in SPM.mat file
%-----------------------------------------------------------------------
try
    save SPMd_MS MS -append
    save SPMd_SS SS -append
catch
    save SPMd_MS MS
    save SPMd_SS SS
end

%=======================================================================
%- E N D: Cleanup GUI
%=======================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')              %-#
spm('FigName','Stats: done',Finter);
fprintf('%-40s: %30s\n','Completed',spm('time'))                      %-#
%fprintf('...use the results section for assessment\n\n')             %-#


%==========================================================================
% function D = get_diag_stat(SPM)
%==========================================================================
function D = get_diag_stat(SPM)
% Build up the menu structure of the all diagnostic and exloratory
% statistics.

xX     = SPM.xX;
DesMtx = xX.X;
nScan  = size(DesMtx,1);

% fMRI-specific operations
if isstruct(xX.K)
    nPS = 0; % Number of scans in previous sessions
    nS = size(DesMtx,1);  % Number of scans
    for s = 1:length(xX.K)
        KH  = xX.K(s).X0;
        nSS = size(KH,1); % Number of scans in this session
        nH  = size(KH,2); % Number of HP bases
        
        %- Zero pad to account for other sessions
        %------------------------------------------------------------------
        KH = [zeros(nPS,nH); full(KH); zeros(nScan-nPS-nSS,nH)];
        DesMtx = [DesMtx KH];
        nPS = nPS + nSS;
    end
end

%- De-mean and remove columns if we don't have a full rank matrix
%--------------------------------------------------------------------------
DWX = DesMtx;
DWX = spm_detrend(DWX);  % Orthogonalize w.r.t intercept
[u, s] = svd(DWX);
DWX = u*s(:,diag(s) > 1e-6);

%- Last ditch sanity check to make sure DWX has the same space as xX.X
%--------------------------------------------------------------------------
DesMtxA = DesMtx;
DesMtxB = [ones(nScan,1) DWX]; %full rank design matrix with intercept

BontoA  = DesMtxA*pinv(DesMtxA)*DesMtxB;

if any(any(abs(DesMtxB-BontoA) > 1e-6))
    disp('Space of Design Matrix Not Matched')
end

%-Load the estimate for global intensity
%--------------------------------------------------------------------------
xGX     = SPM.xGX;
if isempty(xGX.rg)
    VY = SPM.xY.VY;  % Assume handles are valid
    
    %-Compute Global variate
    %=======================================================================
    g = zeros(nScan,1);
    fprintf('%-40s: %30s','Calculating globals',' ')                     %-#
    for i = 1:nScan
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),sprintf('%4d/%-4d',i,nScan)) %-#
        g(i) = spm_global(VY(i));
    end
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')               %-#
else
    g     = xGX.rg;
end
Glob    = [ones(nScan,1) g];


%-Specify all the spatial diagnostic statistics
%-----------------------------------------------------------------------
D(1).name     = 'Corr';
D(1).desp     = 'Durbin-Watson statistic';
D(1).Cr       = [];
D(1).Pdesp    = '-log10(P-value) of Durbin-Watson statistic';
D(1).P        = [];
D(1).diag     = 'dw';
D(1).X        = DWX;
D(1).para     = [];
D(1).help     = ['Assess independence assumption using Durbin-Watson' ...
    ' statistics'];

D(2).name     = 'Norm';
D(2).desp     = 'Shapiro-Wilk statistic';
D(2).Cr       = [];
D(2).Pdesp    = '-log10(P-value) of Shapiro-Wilk statistic';
D(2).P        = [];
D(2).diag     = 'sw';
D(2).X        = [];
D(2).para     = [];
D(2).help     = ['Assess normality assumption using Shapiro-Wilk'...
    ' statistics'];

D(3).name     = 'Homo1';
D(3).desp     = 'Cook-Weisberg Score statistic';
D(3).Cr       = [];
D(3).Pdesp    = '-log10(P-value) of Cook-Weisberg Score statistic';
D(3).P        = [];
D(3).diag     = 'score';
D(3).X        = DesMtx;
D(3).para     = [];
D(3).help     = ['Assess homoscadasticity assumption using Cook-Weisberg' ...
    ' score statistics: With respect to experimental design'];

D(4).name     = 'Homo2';
D(4).Pdesp     = 'Cook-Weisberg Score statistic';
D(4).Cr       = [];
D(4).Pdesp    = '-log10(P-value) of Cook-Weisberg Score statistic';
D(4).P        = [];
D(4).diag     = 'score';
D(4).X        = [];
D(4).para     = 1;
D(4).help     = ['Assess homoscadasticity assumption using Cook-Weisberg' ...
    ' score statistics: With respect to predicted response'];

D(5).name     = 'Homo3';
D(5).desp     = 'Cook-Weisberg Score statistic';
D(5).Cr       = [];
D(5).Pdesp    = '-log10(P-value) of Cook-Weisberg Score statistic';
D(5).P        = [];
D(5).diag     = 'score';
D(5).X        = Glob;
D(5).para     = [];
D(5).help     = ['Assess homoscadasticity assumption using Cook-Weisberg' ...
    ' score statistics: With respect to global signal'];

D(6).name     = 'Dep';
D(6).desp     = 'Cumulative periodogram statistic';
D(6).Cr       = [];
D(6).Pdesp    = '-log10(P-value) of Cumulative periodogram statistic';
D(6).P        = [];
D(6).diag     = 'cp';
D(6).X        = DesMtx;
D(6).para     = [];
D(6).help     = ['Examine long term autocorrelation using cumulative' ...
    ' periodogram statistic'];

D(7).name     = 'Outl';
D(7).desp     = 'Spatial outlier counts';
D(7).Cr       = [];
D(7).Pdesp    = '-log10(P-value) of Spatial outlier counts';
D(7).P        = [];
D(7).diag     = 'outl';
D(7).X        = DesMtx;
D(7).para     = [];
D(7).help     = 'Examine the number of outlier at each voxel';

% D(8).name     = 'ResRMS';
% D(8).desp     = 'Standard deviation of residuals';
% D(8).Cr       = [];
% D(8).Pdesp    = 'Standard deviation of residuals';
% D(8).diag     = 'sdr';
% D(8).X        = DesMtx;
% D(8).para     = [];
% D(8).help     = 'Standard deviation of residuals';
