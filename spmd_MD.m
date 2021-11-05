% The first version is done by Wenlin
% $Id: spmd_MD.m,v 1.36 2008/02/20 12:24:47 nichols Exp $

function varargout = spmd_MD(varargin)
% Create spatial and temporal image around a specific voxel 
% FORMAT spmd_MD(varargin)
% varargin - None
% Using the current spatial and temporal location, show time series plots
% and diagnostic plots for residuals
% _______________________________________________________________________
%
% In the general linear regression, the main tool for the diagnosis of the
% model is to look at the behavior of the residuals. Usually, it is shown
% by graphical inspection of residuals to check the model adequacy. Some
% important diagnostic plots are used here to check the model fitting at
% single voxel. By using this tool, we can find the possible cause of the
% significance of the summary statistics and possible problematic scans. 
% _______________________________________________________________________
% 
% The output of this function includes eight plots. They are:
% Diagnostic plots:
% o Residual vs experimental predictor
% o Residual vs global signal
% o Absolute value of residual vs. experimental predictor
% o Absolute value of residual vs. global signal
% o Lag-1 serial plot
% o Normality plot
% Time series plots:
% o Raw data, drift, and fitted data
% o Standard residual vs horizontal line
% _______________________________________________________________________
%
% Caution:
% Before using this function, there should at least one spatial window
% opened, so that the spatial position could be specified. If bothe
% spatial window exist, this function will automatically select the
% spatial position in the most recently clicked spatial window.
% _______________________________________________________________________
%
% Reference:
% Luo, W-L and Nichols T. E. (2002) Diagnosis and Exploration of
% Massively Univariate fMRI Models. NeuroImage,19:1014-1032, 2003
%
% ______________________________Function Called __________________________
%
%      spm
%      spm_input
%      spm_platform
%      spmd_orthviews
%      spmd_mtsview
%      spmd_getTS
%      spm_figure
%      spmd_MD_plot
%      spm_vol
%      spm_sp
%      spm_filter
%      get_data  (internal)
%      SetWindow (internal)
%_________________________________________________________________________

spm('Pointer','watch');
global TimeCurs


%=============================================================
%-Find the spatial position
%=============================================================
if nargin == 0,
    
    %-Check if both spatial viewers exist and return the handlers
    %---------------------------------------------------------------------
    ms = findobj('Tag','Model Summary');
    sd = findobj('Tag','Scan Detail');

    if isempty(ms) & isempty(sd),
        error('No spatial viewer window exist!');
    end

    %-Select the most recent spatial window, either spatial summary or
    % spatial detail window.
    %---------------------------------------------------------------------
    
    fg = gcbf; 
    
    if isempty(fg) | fg == 1,
        if ~isempty(ms) & ~isempty(sd),
            h   = get(0,'children');
	    tmp = min(find(h==ms),find(h==sd));
	    
            if ~isempty(tmp)
                fg = h(tmp);
            else
                % Shouldn't need this... just in case...
                str = 'Which spatial viewer is the model detail drawn from?';
                input = spm_input(str,1,'bd','Model Summary|Scan Detail',[1,0]);
                if input
                    fg = findobj('Tag','Model Summary');
                else
                    fg = findobj('Tag','Scan Detail');
                end
                close(gcf)
            end
        else
            fg = [ms sd];
        end
    end

    %-get the position of cromshair in the most recently clicked window.
    %--------------------------------------------------------------------
    xyz    = spmd_orthviews('fig',fg,'pos',1);
    xyzmm  = spmd_orthviews('fig',fg,'pos');
    
    if isempty(xyz)
       error('No spatial diagnosis!')
    end
else
    xyz = varargin{1};
end

%- set the toolbar at the top of the model detail window.
%----------------------------------------------------------------------
WS   = spm('WinScale');				%-Window scaling factors
FS   = spm('fontsizes');			% uicontrol font size
PF   = spm_platform('fonts');

%=============================================================
%-Find the temporal position and necessary data
%=============================================================
%-get the temporal position in the temporal summary window. If not
% available, it uses the middle of the scans.
%-------------------------------------------------------------

load SPMd_SS
xSPM  = load(fullfile(SS.SPM.swd,'SPM.mat'));
xX    = SS.xX;
xSPM  = xSPM.SPM;
nScan = size(xX.X,1); 
t     = TimeCurs;
if isempty(t)
    t = nScan/2;
end

%- get all necessary data
%-------------------------------------------------------------
D = get_data(xSPM,SS,xyz);

%=============================================================
%-Fire up the model detail window
%=============================================================
SetWindow(t,xyz,xyzmm);
f    = spm_figure('GetWin','SPMd_MD');

if ~isempty(get(f,'userdata'))

  %- MD window is apparenlty up and running; just re-plot
  %------------------------------------------------------------------

  hs = findobj(f,'Tag','UI_in_MD_plot');
  ht = findobj(f,'Tag','TS_in_MD_plot');
  hit = findobj(f,'Tag','UiTS');
  
  hs = sort(hs);
  ht = sort(ht);
  
  for i=1:length(hs)
    spmd_MD_plot(get(hs(i),'userdata'),get(hs(i),'value'),D);
  end
  
  %- Update time series plots of raw data, drift, fitted data, and
  %  residuals
  %------------------------------------------------------------------
  
  if(get(hit,'value')==1)
      spmd_mtsview('update',ht,... 
                   [D(4:6).data],{D(4:6).desp},...
	        [D(9).data D(3).data],{D(9).desp D(3).desp});
  else
      spmd_mtsview('update',ht,...
		   [D(11:13).data],{D(11:13).desp},...
	       [D(9).data D(10).data ],{D(9).desp D(10).desp});
  end

else  

  %- Create new temporal diagnotic and time series plots
  %------------------------------------------------------------------
  %plot the diagnostic plots

  str = ['Res vs Predictor|Res vs Global Signal|abs(Res) vs Predictor|'...
	 'abs(Res) vs Global Signal|Lag-1 plot|Normality plot'];

  h1 = axes('pos',[0.08 0.745 0.40 0.2],'Visible','on','nextplot','add',...
	    'userdata',D,'Tag','in_MD_plot');
  set(h1,'box','on');
  Ui1 = uicontrol(f,...
		  'style','popup',...
		  'string',str,...
		  'unit','normalized',...
		  'position',[0.08,0.975,0.40,0.02],...
		  'Callback','h = get(gcbo,''userdata'');type = get(gcbo,''value'');spmd_MD_plot(h,type);',...
		  'createfcn','h = get(gcbo,''userdata'');spmd_MD_plot(h,1)',...
		  'value',1,...
		  'FontName',PF.times,'FontWeight','Normal','FontAngle','Normal',...
		  'FontSize',FS(9),...
		  'Interruptible','on','Enable','on',...
		  'ForegroundColor','b',...
		  'userdata',h1,...
		  'Tag','UI_in_MD_plot');

  h2 = axes('pos',[0.55 0.745 0.40 0.20],'Visible','on','nextplot','add', ...
	    'userdata',D,'Tag','in_MD_plot');
  set(h2,'box','on');
  Ui2 = uicontrol(f,'style','popup',...
		  'string',str,...
		  'unit','normalized',...
		  'position',[0.55,0.975,0.40,0.02],...
		  'Callback','h = get(gcbo,''userdata'');type = get(gcbo,''value'');spmd_MD_plot(h,type);',...
		  'createfcn','h = get(gcbo,''userdata'');spmd_MD_plot(h,2)',...
		  'value',2,...
		  'FontName',PF.times,'FontWeight','Normal','FontAngle','Normal',...
		  'FontSize',FS(9),...
		  'Interruptible','on','Enable','on',...
		  'ForegroundColor','b',...
		  'userdata',h2,...
		  'Tag','UI_in_MD_plot');

  h3 = axes('pos',[0.08 0.44 0.40 0.20],'Visible','on','nextplot','add', ...
	    'userdata',D,'Tag','in_MD_plot');
  set(h3,'box','on');
  Ui3 = uicontrol(f,'style','popup',...
		  'string',str,...
		  'unit','normalized',...
		  'position',[0.08,0.67,0.40,0.02],...
		  'Callback','h = get(gcbo,''userdata'');type = get(gcbo,''value'');spmd_MD_plot(h,type);',...
		  'createfcn','h = get(gcbo,''userdata'');spmd_MD_plot(h,5)',...
		  'value',5,...
		  'FontName',PF.times,'FontWeight','Normal','FontAngle','Normal',...
		  'FontSize',FS(9),...
		  'Interruptible','on','Enable','on',...
		  'ForegroundColor','b',...
		  'userdata',h3,...
		  'Tag','UI_in_MD_plot');

  h4 = axes('pos',[0.55 0.44 0.40 0.20],'Visible','on','nextplot','add', ...
	    'userdata',D,'Tag','in_MD_plot');
  set(h4,'box','on');
  Ui4 = uicontrol(f,'style','popup',...
		  'string',str,...
		  'unit','normalized',...
		  'position',[0.55,0.67,0.40,0.02],...
		  'Callback','h = get(gcbo,''userdata'');type = get(gcbo,''value'');spmd_MD_plot(h,type);',...
		  'createfcn','h = get(gcbo,''userdata'');spmd_MD_plot(h,6)',...
		  'value',6,...
		  'FontName',PF.times,'FontWeight','Normal','FontAngle','Normal',...
		  'FontSize',FS(9),...
		  'Interruptible','on','Enable','on',...
		  'ForegroundColor','b',...
		  'userdata',h4,...
		  'Tag','UI_in_MD_plot');
  
  %- time series plots of raw data, drift, fitted data, and residuals
  %- the default is for whitened data
  %- the alternative one is for raw data
  %----------------------------------------------------------------------
  
  
  %- specify the option
  if isfield(xSPM.xVi,'iid') || ...
	(isfield(xSPM.xVi,'form') && strcmp(xSPM.xVi.form,'i.i.d')) || ...
	all(all(xSPM.xVi.V==spdiags(ones(nScan,1),0,nScan,nScan)))
     %- for group data
     str1 = ['Raw       residual'];
     h    = spmd_mtsview('FigYsp',0.6,[],...
                 [D(11:13).data],{D(11:13).desp},...
	       [D(9).data D(10).data ],{D(9).desp D(10).desp});
  else
    %- for sigle subject
    str1 = ['Whitened       residual|Raw       residual'];
    h = spmd_mtsview('FigYsp',0.6,[],...
		   [D(4:6).data],{D(4:6).desp},...
	       [D(9).data D(3).data ],{D(9).desp D(3).desp});
  end
  
  h = h(1:2);
  set(h,'Tag','TS_in_MD_plot');
  set(f,'userdata',[h1 h2 h3 h4]);
  
  
  UiTS = uicontrol(f,'style','popup',...
          'string', str1,...
          'unit', 'normalized',...
          'position',[0.08, 0.37,0.87,0.02],...
          'Callback','h=get(gcbo,''userdata'');type=get(gcbo,''value'');spmd_MDTS_plot(h,type);',...
          'value',1,...
          'userdata',f,...
          'FontName',PF.times,'FontWeight','Normal','FontAngle','Normal',...
          'FontSize',FS(9),...
          'Interruptible','on','Enable','on',...
          'ForegroundColor','b',...
          'Tag','UiTS');
                    
          
  end

spm('pointer','arrow');

%=============================================================
%                 S U B F U N C T I O N S
%=============================================================
function D = get_data(xSPM,SS,xyz)
%-------------------------------------------------------------------------
% Calculate the necessary time series quantities for model detail window. 
%-------------------------------------------------------------------------
global diagV;
global bV;
global vV;
global D;

%- set the variables
%-------------------------------------------------------------------
xX      = SS.xX;
VY      = SS.VY;
Exp     = SS.Exp;
xGX     = xSPM.xGX;
Vbeta   = xSPM.Vbeta;
nScan   = length(VY);
VResMS  = xSPM.VResMS;

%-Check and set data name matrix.
%--------------------------------------------------------------------
if isempty(diagV)
  diagV = VY; % Assume the filehandles are valid 
end
if isempty(bV)
  %-Add full name to Vbeta file
  for i = 1:size(Vbeta,2)
    if isempty(fileparts(deblank(Vbeta(i).fname)))
      Vbeta(i).fname = fullfile(xSPM.swd,Vbeta(i).fname);
    end
  end
  bV = spm_vol(str2mat(Vbeta.fname));
end
if isempty(vV)
  %- Adding full name to VResMS file
  if isempty(fileparts(deblank(VResMS.fname)))
    VResMS.fname = fullfile(xSPM.swd,VResMS.fname);
  end
  vV = spm_vol(VResMS); %-get the variance from ResMS.img
end

%------------------------------------------------------------------
%-Get the data at the current voxel
%-Get all the necessary information for the temporal detail plots,
% including the raw data, drift, residuals,... etc.
%------------------------------------------------------------------
Y    = spmd_getTS(diagV,xyz(1),xyz(2),xyz(3));  % Don't apply gSF! Assume
                                                % already set
beta = spmd_getTS(bV,xyz(1),xyz(2),xyz(3));
sigs = spmd_getTS(vV,xyz(1),xyz(2),xyz(3));

%-Get the high-pass filters matrix and design matrix.
%----------------------------------------------------------------------
X1 = xX.X;  % Just signal of interest
X  = xX.X;  % Full X, to be appended with drift terms
if isstruct(xX.K)
    nPS = 0;                    % Number of scans in previous sessions
    for s=1:length(xX.K)
        KH  = xX.K(s).X0;
        nSS = size(KH,1);       % Number of scans in this session
        nH  = size(KH,2);       % Number of HP bases

        %- Zero pad to account for other sessions
	%-------------------------------------------------------------
        KH     = [zeros(nPS,nH); full(KH); zeros(nScan-nPS-nSS,nH)];
        X = [X KH];
        nPS    = nPS + nSS;
    end
end

X0 = X(:,size(X1,2)+1:end);  % Just drift terms

%----------------------------------------------------------------------
%
% This should work for drift terms
%
%  Y = [X1  X0][b1' b0' ]' + \eps
%    = [X10 X0][b1' b0*']' + \eps            X10 is X1 orthog wrt X0
%
% This is really what happens in spm_spm:
%
%    Res = KWY      - KWX10 b1
%        = Y - X0 b0* - X10 b1
%        = Y - X0 b0  - X1  b1
%
% Hence, to get the fitted values associated with the drift...
%
%        X0 b0  = Y - Res - X1 b1
%
%
% Doesn't work for coloring/temporal smoothing!
%----------------------------------------------------------------------
W     = xX.W;                    %-Get whitening/Weighting matrix

H2    = W*X*pinv(W*X); %- Calculate the Hat matrix for whitened data
H3    = X*pinv(W*X)*W; %- Calculate the Hat matrix for raw data

I     = speye(length(W));
ResW  = (I-H2)*W*Y;              %- Whitened residual
Res   = (I-H3)*Y;                %- Raw residual

h2ii  = diag(H2);                %- Diagonal matrix from the Hat matrix for whitened data
h3ii  = diag(H3);                %- Diagonal matrix from the Hat matrix for the raw data


sResW  = ResW./sqrt(sigs*(1-h2ii)); %- Studentized standard error for whitened residual
sRes   = Res./sqrt(sigs*(1-h3ii));  %- Studentized standard error for raw residual

I2 = (1-h2ii)<sqrt(eps);   %-Catch problems with 'negative' variance 
I3 = (1-h3ii)<sqrt(eps);   % due to dummy variables negating bad observations
sResW(I2) = 0;
sRes(I3) = 0;

KWY   = spm_filter(xX.K,W*Y);      % Whitened, De-drift-ed data for whitened data
KY    = spm_filter(xX.K,Y);        % De-drift-ed data for raw data

KWX   = spm_sp('Set',spm_filter(xX.K,W*X1)); %- Structure of K*W*X for whitened data
KX    = spm_sp('Set',spm_filter(xX.K,X1));   %- Structure of K*X for raw data

%- raw data
Yh1  = X1*beta;                    % Experimental fit only
Yh0  = Y - Res - Yh1;              % Drift fit only
Yh   = Yh1 + Yh0;                  % Fit
Yh0  = Yh0 + mean(Yh1);            % Move mean to Yh2
Yh1  = Yh1 - mean(Yh1);

%- whitened data
WY   = W*Y;                        
Yh1W = W*X1*beta;                  % Experimental fit only
Yh0W = W*Y - ResW - Yh1W;          % Drift fit only
YhW  = Yh1W + Yh0W;                % Fit
Yh0W = Yh0W + mean(Yh1W);          % Move mean to Yh2W
Yh1  = Yh1W - mean(Yh1W);         


LRes = nScan;

%----------------------------------------------------------------------
%set up the data structure for the residual plots and time series plots.
%----------------------------------------------------------------------
% D.data = Vector of data
% D.desp = Description of data
% D.help = Tooltip help string
% D.abline = abline specification; e.g. [0 1] plots the identity line
%            'h=0' gives horizontal line at zer0, 'v=2' gives vertical
%            line at 2, 'hv=0' gives horizontal and vertical lines.

D(1).data = ResW;
D(1).desp = 'Residual';
D(1).help = 'Residuals calculated from the specified model (including possible whitening and de-drifting)';
D(1).abline= 'h=0';  

D(2).data = abs(ResW);
D(2).desp = '|Residual|';
D(2).help = 'Absolute value of residuals calculated from the specified model (including possible whitening and de-drifting)';
D(2).abline= [];  

D(3).data = sResW;
D(3).desp = 'Std. Residual (W)';
D(3).help = 'Studentized Residuals from the specified model (including possible whitening and de-drifting)';
D(3).abline= 'h=0';  

D(4).data = WY;
D(4).desp = 'Data';
D(4).help = 'Raw data after whitened';
D(4).abline= []; 

D(5).data = Yh0W;
D(5).desp = 'Drift';
D(5).help = 'Fitted drift in the fMRI time series';
D(5).abline= [];  

D(6).data = YhW;
D(6).desp = 'Fitted';
D(6).help = 'Fitted data of full model, including the drift.';
D(6).abline= [];  

D(7).data = Exp.Pred;
D(7).desp = 'Predictor';
D(7).help = 'Experimental Predictor in the model';
D(7).abline= [];  

% D(8).data = xGX.rg;
D(8).data = SS.GX.ts; % Tom's change
D(8).desp = 'Global Signal';
D(8).help = 'Global signal over the space';
D(8).abline= [];  

D(9).data = zeros(1,nScan)';
D(9).desp = 'HZ Line'      ;
D(9).help = '';
D(9).abline= [NaN NaN];  

D(10).data = sRes;
D(10).desp = 'Std. Residual';
D(10).help = 'Studentized Residuals from the specified model (including possible whitening and de-drifting)';
D(10).abline= 'h=0';  

D(11).data = Y;
D(11).desp = 'Data';
D(11).help = 'Raw data after whitened';
D(11).abline= []; 

D(12).data = Yh0;
D(12).desp = 'Drift';
D(12).help = 'Fitted drift in the fMRI time series';
D(12).abline= [];  

D(13).data = Yh;
D(13).desp = 'Fitted';
D(13).help = 'Fitted data of full model, including the drift.';
D(13).abline= [];  

return

function h = SetWindow(t,xyz,xyzmm)
%----------------------------------------------------------------------
% SetWindow: set up the temporal detail window.
%
% FORMAT SetWindow(t,xyz)
%  t: current temporal position
%  xyz: most recent spatial position
%  The output is temporal detail window, Tagged "SpatempFig"
%----------------------------------------------------------------------
h = findobj('Tag','SPMd_MD');
if isempty(h)
    h = figure(...
        'Tag','SPMd_MD',...
        'numbertitle','off',...
        'unit','normalized',...
        'Position',[0.1 0 0.5 0.85],...
        'MenuBar','none',...
        'windowbuttondownfcn','h = get(gcf,''userdata'');t = spmd_pointer(h);spmd_mtsview(''set'',t);');
else
    figure(h)
end
set(h,'Name',sprintf('Model Detail: Voxel %d,%d,%d (%.1fmm,%.1fmm,%.1fmm)  Scan %d',...
		     round(xyz),xyzmm,round(t)));

return

