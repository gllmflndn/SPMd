function spmd_SS(varargin)
% Create page for scan summary plots
% FORMAT spmd_SS
% By default, shows the following plots
%      Experimental Predictor
%      Global intensity
%      Spatial Outlier count
%      Registration shift movements (if available)
%      Registration rotation movements (if available)
%      Periodogram
%
% FORMAT spmd_SS(Y1,N1,Y2,N2,...)
% Specify each timeseries manually
% Y*  - n-by-p matrix
% N*  - p-length cell array containing the name of each column of Y1
%
% FORMAT spmd_SS('Reposition',tr)
% Move the temporal cursor (only for existing plot)
% tr  - scan to be focused on
%
% FORMAT spmd_SS('TR',tr)
% Set initial position of temporal cursor to be tr (must be first
% argument).
% tr  - scan to be focused on
%__________________________________________________________________________
%
% This is essentially a wrapper for spmd_mtsview
%__________________________________________________________________________
% @(#)spmd_SS.m 1.19 Tom Nichols & Wen-Lin Luo 03/07/16


%-open the temporal summary window
SetWindow;

%-Add print and close button at the top of the window
%--------------------------------------------------------------------------

Marker     = 'o';
MarkerSize = 3;

fg = spm_figure('GetWin','SPMd_SS');
FS = spm('fontsizes');

uicontrol(fg,...
    'String','Print',...
    'Units','normalized',...
    'Position',[0.01 0.97 0.15 0.03],...
    'CallBack','spm_figure(''Print'',gcbf)',...
    'FontSize',FS(10),...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Tag','ToolBar',...
    'HandleVisibility','callback',...
    'ForegroundColor','k');

uicontrol(fg,...
    'String','Close',...
    'Units','normalized',...
    'Position',[0.17 0.97 0.15 0.03],...
    'CallBack','close(gcbf)',...
    'FontSize',FS(10),...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Tag','ToolBar',...
    'HandleVisibility','callback',...
    'ForegroundColor','k');

%==========================================================================
%-Specify the position of temporal cursor
%==========================================================================
TR = [];

if ~nargin
    fNm = {'predint','global','toutlier','shift','rotate','pg'};
end

for i=1:2:numel(varargin)
    switch lower(varargin{i})
        
        case 'ts'
            fNm = varargin{i+1};
            
        case 'tr'
            TR  = varargin{i+1};
            
        case {'reposition', 'repos'}
            spmd_mtsview('Set',varargin{i+1});
            return;
    end
end

%==========================================================================
%-Assign the data
%==========================================================================

D = diag_SS;
for j = 1:numel(fNm)
    Id(j) = strmatch(fNm(j),{D(:).fname});
end

Id = sort(Id);

TS = {};
for i = 1:length(Id)
    if Id(i) ~= 6
        TS  = {TS{:} D(Id(i)).ydata {D(Id(i)).title D(Id(i)).legend}};
    end
end

nTS = length(Id);

%==========================================================================
%-Create the time series plots
%==========================================================================
if any(Id == 6)
    if nTS > 1
        spmd_mtsview('FigYsp',[0 (nTS-1)/nTS],TR,TS{:});
    end
    AX = axes(...
        'Position',[0.08 0.05 0.87 1/nTS-0.07],...
        'Visible','off',...
        'FontSize',FS(8));
    plot(AX,D(6).xdata, D(6).ydata ,'Marker',Marker,'MarkerSize',MarkerSize);
    title(AX,D(6).title,'fontsize',FS(8));
    xlabel(AX,'Frequency (Hz)','FontSize',FS(8));
else
    spmd_mtsview(TR,TS{:});
end


function SetWindow
%==========================================================================
%-Create the temporal summary window.
%==========================================================================

h = findobj('Tag','SPMd_SS');

WS   = spm('WinScale');

Rect = spm('WinSize','Graphics','raw').*WS; %-Graphics window rectangle
Rect = round(Rect);
S    = get(0,'screensize');
FW   = Rect(4)/2;
FL   = 2*Rect(4)/3;

if isempty(h)
    figure(...
        'Tag','SPMd_SS',...
        'Units','normalized',...
        'Position',[0.1 0.27 0.5 0.7],...
        'MenuBar','none',...
        'NumberTitle','off',...
        'Name',sprintf('Scan Summary: Scan ?'));
else
    figure(h)
end


function D = diag_SS(varargin)
%==========================================================================
%-Preparation for global signal window
%-include global signal, drift, and drift+covariates in the global
% signal window.
%==========================================================================

load SPMd_SS;
load SPM;
PG     = SS.PG;
GX     = SS.GX;
Toutl  = SS.Toutl;
Exp    = SS.Exp;
RX     = SS.RX;
%RT    = SPM.xY.RT;
%Tom's change
if isfield(SPM.xY,'RT')
    RT = SPM.xY.RT;
else
    RT = 1;
end
%%%%%%%%%%%%%%%%%%%
nScan  = length(SS.VY);
Xdata  = (1:nScan)';


if isfield(SS,'Exp')
    D(1).ydata = SS.Exp.Pred;
    D(1).xdata = Xdata;
else
    D(1).data  = [];
end
D(1).title = 'Experimental Predictor';
try
    D(1).legend = SS.Exp.PredNms;
catch
    D(1).legend = {};
end

D(1).fname = 'predint';
D(1).help  = 'The predicted signal';

str = sprintf('Global Signal: Ftest=%0.2f, P-value=%0.2g',GX.Fstat,GX.P);

D(2).ydata  = [GX.ts GX.Est];
D(2).xdata  = Xdata;
D(2).title  = str;
D(2).legend = {'Global' 'Drift' 'Drift+Fit'};
D(2).fname  = 'global';
D(2).help   = 'The global intensity';

D(3).ydata  = Toutl.prop;  %- change to proportion
D(3).xdata  = Xdata;
D(3).title  = 'Spatial Outlier Rate';
D(3).legend = {};
D(3).fname  = 'toutlier';
D(3).help   = 'Percent of expected outliers relative to nominal \alpha=0.05';

if isfield(SS,'RX') && ~isempty(SS.RX)
    D(4).ydata = [RX(1).ts RX(2).ts RX(3).ts];
    str1       = {sprintf('Shift Movement: F=(%0.2f,%0.2f,%0.2f), P=(%0.2f,%0.2f,%0.2f)',...
        RX(1).Fstat,RX(2).Fstat,RX(3).Fstat,RX(1).P,RX(2).P,RX(3).P)};
    D(5).ydata = [RX(4).ts RX(5).ts RX(6).ts];
    str2       = {sprintf('Rotation Movement: F=(%0.2f,%0.2f,%0.2f), P=(%0.2f,%0.2f,%0.2f)',...
        RX(4).Fstat,RX(5).Fstat,RX(6).Fstat,RX(4).P,RX(5).P,RX(6).P)};
else
    D(4).ydata = [];
    str1       = [];
    D(5).ydata = [];
    str2       = [];
end

D(4).xdata  = Xdata;
D(4).title  = str1;
D(4).legend = {'X','Y','Z'};
D(4).fname  = 'shift';
D(4).help   = 'Registration shift movement parameters';

D(5).xdata  = Xdata;
D(5).title  = str2;
D(5).legend = {'pitch','roll','yaw'};
D(5).fname  = 'rotate';
D(5).help   = 'Registration rotation movement parameters';

D(6).xdata  = PG.freq/(2*RT);
D(6).ydata  = PG.power;
D(6).title  = 'Average Periodogram';
D(6).legend = {'Periodogram' 'Frequency'};
D(6).fname  = 'pg';
D(6).help   = 'Average periodogram of raw residuals';
