% The first version is done by Wenlin
% $Id: spmd_comp_SS.m,v 1.22 2008/02/20 12:22:54 nichols Exp $

function spmd_comp_SS(varargin)
% Compute all necessary variables for scan summary
% FORMAT spmd_comp_SS(varargin)
% varargin - names of variables that users want to calculate
% _______________________________________________________________________
%
% The output of variables are appended to SPMd_SS
% 
% Predictor                   - indicator of predicators of interest
% Global Signal               - matrix of global signal related time
%                               series. Three colums with length equal to
%                               the number of scans in this matrix:
%                               global signal, predicted global signal by
%                               all experimental covariates, and
%                               predicted global signal by nuisance covariates.
% Prop. of Toutliers          - The proportion of the number of outlier
%                               in the actual experiment to the expected outlier
% Reg. Parameters             - matrix of registration parameters
% Average periodogram         - average periodogram of raw residuals
% _______________________________________________________________________
%
% Warnings
% 1. The realignment file must be there.
% 2. Only "ONE" file for realignment parameters is allowed here. For
%    multiple sessions or multiple subjects experiments, please
%    concatenate them first!
% 3. Before using this function, users need to run SPM and spmd_comp_MS.m
%    first. 
% _______________________________________________________________________
%
% Reference:
% Luo, W-L and Nichols T. E. (2002) Diagnosis and Exploration of
% Massively Univariate fMRI Models. NeuroImage,19:1014-1032, 2003
%
% ______________________________Function Called __________________________
%      spm_str_manip
%      spm_select
%      spm_Ncdf
%      spm_vol
%      spm_global
%      spm
%      spm_input
%      spm_Fcdf
%      spmd_comp_MS
%      GetExpPred (internal)
%      GetTOutl   (internal)
%      GetGlob    (internal)
%      GetPG      (internal)
%      GetRegParm (internal)
%      GetF       (internal)
%_________________________________________________________________________



spm_defaults
%-Check if SPMd_SS.mat in current directory
%------------------------------------------------------------------------
if exist(fullfile('.','SPMd_SS.mat'),'file')==0
  swd = '.';
  if exist(fullfile(swd,'SPM.mat'))~=2,
    error('No SPM.mat in current directory')
  end
  xSPM = load(fullfile(swd,'SPM.mat'));
  xSPM = xSPM.SPM;
  xX   = xSPM.xX;
  VY   = xSPM.xY.VY;
  if isfield(xSPM,'Sess'), Sess = xSPM.Sess; else, Sess = []; end

  [nScan nVar]  = size(VY);    
  DesMtx        = xX.X;

  %- fMRI-specific operations
  %----------------------------------------------------------------------
  if isstruct(xX.K)
    
    nPS = 0;              % Number of scans in previous sessions
    nS  = size(DesMtx,1); % Number of scans
    for s=1:length(xX.K)
      KH = xX.K(s).X0;
      nSS = size(KH,1);   % Number of scans in this session
      nH  = size(KH,2);   % Number of HP bases
                          % Zero pad to account for other sessions
      KH = [zeros(nPS,nH); full(KH); zeros(nScan-nPS-nSS,nH)];
      DesMtx = [DesMtx KH];
      nPS = nPS + nSS;
    end
  end
  S   = xSPM.xVol.S;       %- number of voxels in each scan
  fDM = DesMtx;
  rDM = xX.X;
   
else
  load SPMd_SS
  swd = SS.SPM.swd;
  if exist(fullfile(swd,'SPM.mat'),'file')~=2 
     spmmat_handle = spm_select(1,'SPM.mat','Select SPM.mat file');
     swdSPM = spm_str_manip(spmmat_handle,'hd');
     xSPM = load(fullfile(swdSPM,'SPM.mat'));
 else
     xSPM = load(fullfile(swd,'SPM.mat'));
 end
  
  xSPM = xSPM.SPM;
  xX  = SS.xX;
  VY  = SS.VY;
  if isfield(xSPM,'Sess'), Sess = xSPM.Sess; else, Sess = []; end
  [nScan nVar]  = size(VY);   
  S   = xSPM.xVol.S;       %- number of voxels in each scan
  rDM = xX.X;
  fDM = xX.fDM;
  swd = pwd;
end

%-If there is no high pass filter, then the reduced model is model with
% only intercept term.
%------------------------------------------------------------------------
if size(rDM,2) == size(fDM,2),
  rDM = ones(nScan,1);
end

xX.fDM  = fDM;
xX.rDM  = rDM; 

SS.xX   = xX;
SS.VY   = VY;
SS.RX   = [];

%- Action to take
%------------------------------------------------------------------------

if nargin == 0,
  varargin = {'predint','toutlier','global','pg','regparm'};
end

a = 1;

while (a <= length(varargin))
  switch lower(varargin{a})

   case 'predint'
    [ExpPred,PredInt,PredNms] = GetExpPred(xSPM,[]); % predictors of interest
    Exp.name    = 'Predictors of interest';
    Exp.Ipi     = PredInt;
    Exp.Pred    = ExpPred;
    Exp.PredNms = PredNms;
    SS.Exp      = Exp;
    a = a+1;
       
   case 'toutlier'
    Toutlier     = GetTOutl(SS,xSPM);                  % Outlier Count
    Toutl.name   = 'Spatial outlier rate (% of nominal)';
    Toutl.ts     = Toutlier;
    Toutl.prop   = (Toutlier/(2*spm_Ncdf(-3)*S))*100;  % The percent of the number of
                                                       % outliers obs of number expected
    SS.Toutl     = Toutl;
    a = a+1;
   
   case 'global'
    glob         = GetGlob(SS);                        % Global signal

    %-Calculate the F-test for the correlation between the temporal
    % measurements to the predictors.
    %------------------------------------------------------------------------------------
    [ts,Stat]	  = GetF(glob,rDM,fDM);
    GX.name       = 'Temporal Global Signal';
    GX.ts         = ts(:,1);
    GX.Est        = ts(:,2:3);
    GX.Fstat      = Stat(1);
    GX.P          = Stat(2);
    SS.GX         = GX;
    
    a = a+1;
    
   case 'pg'
    [power,freq]   = GetPG(SS,xSPM);                   % Periodogram
    PG.name        = 'Periodogram';
    PG.power       = power;
    PG.freq        = freq;
    SS.PG          = PG;
    a = a+1;
   
   case 'regparm'
    RegParm     = GetRegParm(VY);                      % Registration parameters
    
    %-Calculate the F-test for the correlation between the temporal
    % measurements to the predictors.
    %---------------------------------------------------------------------------------------
    if ~isempty(RegParm)
      RegName     = char('x-shift','y-shift','z-shift','pitch','roll','yaw');
      ParmI       = 1:6;
     % close(gcf)
      ParmI       = sort(ParmI);
      for i = 1:length(ParmI),
	[ts,Stat]   = GetF(RegParm(:,ParmI(i)),rDM,fDM);
	RX(i).name  = RegName(ParmI(i),:);
	RX(i).ts    = ts(:,1);
	RX(i).Est   = ts(:,2:3);
	RX(i).Fstat = Stat(1);
	RX(i).P     = Stat(2);
	clear ts Stat
      end
    else
      RX            = [];
    end
    SS.RX           = RX;
    a = a+1;
        
   otherwise
    error('Unknown command');
  end
end

SS.SPM.swd = swd;

try
  save SPMd_SS SS -append
catch
  save SPMd_SS SS
end
%close(gcf)
return

%=======================================================================
%                      S U B F U N C T I O N S
%=======================================================================
function [X1, Ipi, PredNms] = GetExpPred(SPM,Con,ShCov)
%======================================================================
%  Get experimental predictors, indicator of predictors of interest, and
%  names 
%   SPM - SPM structure
%   Con - Contrast of interest
%   ShCov - Show covariates?  (Defaults to 0, no)
%
%  Rules:
%   o If using hrf basis (alone or with derivatives) only show canonical,
%     and don't user-specified covariates.
%   o If using any other basis, give up, and show everything, except
%     user-specified covariates
%   o If *no* conditions and only user-specified covariates, show covariates
%  Note: ShCov flag overrides this, and shows the covariates regardless
%======================================================================
if isfield(SPM,'Sess'), Sess = SPM.Sess; else, Sess = []; end
xX   = SPM.xX;
X    = xX.X;
p    = size(X,2);
try 
  xBFnm = SPM.xBF.name;
catch
  xBF = '';
end
if nargin<2, Con = []; end
if nargin<3, ShCov = 0; end

VisPred  = logical(ones(1,p));          % Columns to show
%- for single subject analysis
if (length(VisPred)>1)
  VisPred(all(diff(X)==0)) = 0;  % Hide any constant preds
end

% Build up predictor names

% Default names
PredNms = cell(1,p);
for i=1:p
  PredNms{i} = strrep(xX.name{i},'Sn(1) ','');
end
% Can we do better?
if ~isempty(Sess)
  if length(Sess.U)>0
    if strncmp(SPM.xBF.name,'hrf',3)
      % Pick out one predictor per condition; don't show derivatives, etc
      for i=1:length(Sess.Fc)
	for j=1:length(Sess.Fc(i).i)
	  if j==1
	    PredNms{Sess.Fc(i).i(j)} = Sess.Fc(i).name; 
	  else
	    VisPred(Sess.Fc(i).i(j)) = 0;   % Hide secondary columns
	  end
	end
      end
    else
      % Don't know how to deal with other basis types... use def nms, all col
    end
    % Hide user-specified covariates
    %if ~ShCov
    %  VisPred(SPM.Sess.col(end-size(Sess.C.C,2):end)) = 0;
    %end
    
    if (~ShCov & size(Sess.C.C,2)>0)
      %VisPred(Sess.col(end-size(Sess.C.C,2):end)) = 0;
      colnum = size(Sess.col,2) - size(Sess.C.C,2);
      VisPred(Sess.col(end-size(Sess.C.C,2) + colnum :end)) = 0;
    end
  else
    % Condition only design; use default names, all cols
  end
else
  % Non-fMRI design; use default names, all cols
end

PredNms = PredNms(VisPred);

if isempty(Con)
  
  % No contrast supplied
  if isempty(VisPred)

    % No non-constant predictors!
    X1 = X;
    Ipi = logical(ones(1,size(X,2)));  % Should just be logical([1])

  else

    %- if there are more than one column,
    %- For each row, find maximum value and set all other elements in this
    %  row to NaN, and make sure every row must have exactly one non_NaN
    %  entry
    
    % X1 = X(:,VisPred);
    % X1(repmat(sf_nanmax(X1')',1,size(X1,2))~=X1)=NaN;
    % X1(X1==0) = NaN; %- Change 0's to NaN;
        
    X1 = X(:,VisPred);
    if size(X1,2)>1
        X1(repmat(sf_nanmax(X1')',1,size(X1,2))~=X1)=NaN;
        X1(X1==0) = NaN; %- Change 0's to NaN;
    end

    % Now deal w/ all NaN rows
    % find first row that isn't all zeros (NaNs now)
    tmp = min(find(~isnan(sf_nanmax(X1')'))); 
    if ~isempty(tmp)
      firstNonNaN = min(find(~isnan(X1(tmp,:))));
      % Now proceed down the rows, choosing the column based on the last
      % previous non-NaN column
      for i = find(all(isnan(X1)'))
	if ~isempty(firstNonNaN)
	  X1(i,firstNonNaN) = 0;
	  firstNonNaN = [];
	else
	  X1(i,min(find(~isnan(X1(i-1,:))))) = 0;
	end
      end
    end
    
    % % Scaling is good for DesMtx images, but not sure if it's so good for
    % % plotting... 
    % X1 = spm_DesMtx('sca',X1(:,VisPred));

    Ipi  = zeros(p,1);
    Ipi(VisPred) = 1;
    Ipi = logical(Ipi);
  end

else

  % Contrast supplied

  % Message Con
  Con = [Con zeros(1,p-length(Con))];

  Ipi  = zeros(p,1);
  Ipi(Con~=0) = 1;
  Ipi = logical(Ipi);
  
  X1 = X*Con';

  PredNms = {'Contrast'}
end

return


function TOutlier = GetTOutl(SS,xSPM)
%======================================================================
%  Get temporal outlier profile
%======================================================================

if ~isfield(SS,'Toutl')
  TOutlier = [];
else
  TOutlier = SS.Toutl.ts;
end
    
if isempty(TOutlier)
  spmd_comp_MS(xSPM,{'Outl'});
  load SPMd_SS
  TOutlier = SS.Toutl.ts;
end
return
   

function Glob = GetGlob(SS)
%======================================================================
%  Get temporal outlier profile
%====================================================================== 
if ~isfield(SS,'GX')
  Glob = [];
else
  Glob = SS.GX.ts;
end

if isempty(Glob)
  VY    = SS.VY; % Assume handles are valid
  q     = length(VY);
  Glob  = zeros(q,1);      

  fprintf('%-40s: %30s','Calculating globals',' ')                     %-#
  for i = 1:q
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),sprintf('%4d/%-4d',i,q)) %-#
    Glob(i) = spm_global(VY(i));
  end
end

return

function [power,freq] = GetPG(SS,xSPM)
%======================================================================
%  Get average periodogram of raw residuals
%======================================================================
if ~isfield(SS,'PG')
  power = [];
  freq  = [];
else
  power = SS.PG.power;
  freq  = SS.PG.freq;  
end
if isempty(power)
  spmd_comp_MS(xSPM,{'Dep'});
  load SPMd_SS
  power = SS.PG.power;
  freq  = SS.PG.freq;  
end

return


function RegPar = GetRegParm(VY)
%======================================================================
%  Get registration parameters if possible
%======================================================================
nScan = size(VY,1);

if (nargin<2 | isempty(ParmfNm)) 
  fNm    = VY(1).fname;
  Path   = spm_str_manip(fNm,'H');
  BaseNm = spm_str_manip(fNm,'tr');
  
  if (BaseNm(1)=='n'), BaseNm(1) = []; end
  if (BaseNm(1)=='r'), BaseNm(1) = []; end   % Get down to raw filename
  
  ParmfNm = [ Path '/rp_' BaseNm '.txt' ];     

end

if ~exist(ParmfNm)
   str = {'Data directory does not contain realignment parameters file:',...
  	 ['rp_' BaseNm '.txt or .dat' ],['(dir = ',Path,')'],' ',...
  	 'Do you want to specify the file?'};
    
  input = spm_input(str,1,'bd','Yes|No',[1,0]);
  if input, 
    ParmfNm = spm_select(1,'.*\.txt$','Select realignment file',Path);
  else
    ParmfNm = [];
    fprintf('%-40s: %30s\n\n',...
	    'No realignment parameter file specified',spm('time')); 
    RegPar = [];
  end
end

if ~isempty(ParmfNm)
  RegPar = load(ParmfNm); 
  
  %-Check if the length of the realignment parameters is the same as the
  %-nScan
  %---------------------------------------------------------------------
  LRegPar = size(RegPar,1);
  if LRegPar == nScan,
    RegPar = RegPar;
  elseif LRegPar < nScan
    error(['The length of realignment parameters is less than the number' ...
	   ' of scan: Wrong file is specified!']);
  else
    str1 = sprintf('Need realignment parameters for %0.0f scans,',nScan);
    str2 = sprintf('but %0.0f rows are found!',LRegPar);
    spm_input({str1,str2},1,'bd','OK');    
    while LRegPar ~= nScan,
      tmp = spm_input('Enter unnecessary scans','+1','n');
      if max(tmp) > LRegPar;
	str3 = ['The indices of omitted scans exceed the length'...
		' of realignment parameter!'];
	spm_input(str3,1,'bd','OK');
      else
	RegPar = load(ParmfNm);
	RegPar(tmp',:)=[];
	LRegPar = size(RegPar,1);
	if LRegPar ~= nScan,
	  str1 = sprintf('Need realignment parameters for %0.0f scans,', ...
			 nScan);
	  str2 = sprintf('but %0.0f rows are found!',LRegPar);
	  spm_input({str1,str2},1,'bd','OK');    
	end
      end
    end
  end
  
  % Nuke first if a reference
  if all(RegPar(1,:)==[0 0 0 0 0 0])
    RegPar(1,:) = NaN;
  end
end
return

function [Global,Stat] = GetF(glob,X1,X)
%======================================================================
%-Calculate the F-test for the correlation between the temporal
% measurements to the predictors
%   X1:   matrix with less variables 
%   X :   matrix with more variables, the first nVar columns contains the
%         variables in X1.
%======================================================================
[nScan nVar1] = size(X1);
nVar          = size(X,2);
X2            = X(:,nVar1:end);
R             = glob-X*pinv(X)*glob;     % residuals for Full model
R2            = glob-X2*pinv(X2)*glob;   % residuals for Reduced model
RSS           = sum(R.^2);               % Sum of Squared residuals for Full model
RSS2          = sum(R2.^2);              % Sum of Squared residuals for Reduced model
%F             = (RSS2-RSS)/((nVar1-1)*RSS)*(nScan-nVar);

% using for contrast
if ((nVar-nVar1)==0),
    F = 1;
    p = 1-spm_Fcdf(F,1,nScan-nVar); 
else
% using for normalsituation
    F = (RSS2-RSS)/((nVar-nVar1)*RSS)*(nScan-nVar);
    p = 1-spm_Fcdf(F,nVar-nVar1,nScan-nVar); 
end
% p             = 1-spm_Fcdf(F,nVar1-1,nScan-nVar); 
Ghni          = glob-R2;
Gh            = glob-R;
Global        = [glob Ghni Gh];
Stat          = [F p];
return

% Sub-function: nanmax
function m = sf_nanmax(A)

% Assume that no column is all NaN's

A(isnan(A)) = -Inf;
m = max(A);
nanIndex = find(m==-Inf);
m(nanIndex) = NaN;
return
