function varargout = spmd_mtsview(TR,varargin)
% Create viewer of multiple time series.
% FORMAT spmd_mtsview('set',t)
% Input:
%  t          - temporal position of the cursor.
% FORMAT h = spmd_mtsview('get')
%  h (output) - current temporal position
% FORMAT spmd_mtsview('clickset')
%  Move the cursor to the clicked temporal position
% FORMAT spmd_mtsview('FigYsp',FigPos,TR,TS1[,TS2,...])
%  FigPos     - Position of the time series plots in the figure window.
%  TR         - Temporal position
%  TS1        - Time series
% FORMAT spmd_mtsview('update',AxHd,TS1[,TS2,...])
%  AxHd       - Handles of axes to update
%  TS1,TS2... - Time series (Note, there must one for each row of AxHd)
% FORMAT h = spmd_mtsview(TR,TS1[,TS2,...])
%  TR         - Temporal position
%  TS1        - Time series
%  h (output) - axes handle of the time series plots.
%__________________________________________________________________________
% To move the temporal cursor to timepoint 13, use the following form
%
%    spmd_mtsview('Set',13)
%__________________________________________________________________________
% %W% Tom Nichols %E%

%----------------------------- Function called --------------------------
%      spm
%      spm_platform
%      spmd_pointer
%      SetTimeCurs (internal)
%      GetTimeCurs (internal)
%      ClickMove (internal)
%      countTS (internal)
%      InitTimeCurs (internal)
%      MakeAxesPos (internal)
%------------------------------------------------------------------------

global SPMd_defs;

FigPos = [];

if (~isempty(TR))
    if isstr(TR)
        switch (lower(TR))
            case 'set'
                SetTimeCurs(varargin{1});
                return
            case 'get'
                varargout{1} = GetTimeCurs;
                return
            case 'clickset'
                ClickMove;
                return
            case 'figysp'
                FigPos = varargin{1};
                TR = varargin{2};
                varargin(1:2) = [];
            case 'update'
                UpdateTS(varargin{1:end});
                return
            otherwise
                error('Unknown command');
        end
    end
end

%-window Parameters
%-----------------------------------------------------------------------
WS     = spm('WinScale');
FS     = spm('FontSizes');
PF     = spm_platform('fonts');

[TSs,Titl,Legs,nTS,nTR] = countTS(varargin{:});

InitTimeCurs(nTR);
Tcurs = GetTimeCurs;

AxPos = MakeAxesPos(nTS,FigPos);

XlimTR  = [0.5 nTR+0.5];
if ~isempty(TR), XlimSec = ([0.5 nTR+0.5]-1)*TR; end

AxHnd = zeros(nTS,2);
MarkerSize = SPMd_defs.MarkerSize;
Marker      = SPMd_defs.Marker;

for i = 1:nTS
    TS = TSs{i};
    
    AxHnd(i,1) = axes('Position',AxPos{i},...
        'box','on',...
        'ButtonDownFcn','spmd_mtsview(''ClickSet'')');
    
    hCur  = line(Tcurs*[1 1],[NaN NaN],...
        'color','red','LineStyle',':',...
        'Tag','TimeCursor',...
        'HitTest','off');
    hLin  = line(1:nTR,TS,'HitTest','off','marker',Marker,'markersize',MarkerSize);
    
    set(hLin,'Tag','Temporal');
    set(gca,'xlim',XlimTR,'UserData',hCur,'fontsize',FS(8));
    
    if (~isempty(Titl{i}))
        ht = title(Titl{i});
        set(ht,'fontsize',FS(8));
    end
    
    if (i==nTS)
        hx = xlabel('time (TRs)');
        set(hx,'fontsize',FS(8));
    else
        set(gca,'XTickLabel',[]);
    end
    
    if ~isempty(Legs{i})
        hl = legend(hLin,Legs{i});
        set(hl,'fontsize',FS(8),'Tag','legend','HitTest','off');
    end
    set(hCur,'Ydata',get(gca,'Ylim'));
    
    if ~isempty(TR)
        %- Add top (seconds) axis
        %----------------------------
        AxHnd(i,2) = axes('Position',get(AxHnd(i,1),'Position'),...
            'Xlim',XlimSec,...
            'Ylim',get(AxHnd(i,1),'Ylim'),...
            'YTickLabel',[],...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');
        
        if (i==1)
            xlabel('time (secs)');
        else
            set(AxHnd(i,2),'XTickLabel',[])
        end
        
    end
    
end

if nargout
    varargout = {AxHnd};
end


function [TS,Titl,Leg,nTS,nTR] = countTS(varargin)
%---------------------------------------------------------------------------
% Deal out arguments into variables.  Skip empty time series and their
% legends.
%---------------------------------------------------------------------------
nArg = length(varargin);

TS   = cell(1,nArg);
Titl = cell(1,nArg);
Leg  = cell(1,nArg);

k = 0;
i = 1;

while (i <= nArg)
    if isempty(varargin{i})
        if (i<nArg) && iscell(varargin{i+1})
            i = i+1;
        end
    else
        
        k = k+1;
        TS{k} = varargin{i};
        
        if (i==1)
            nTR = length(TS{k});
        else
            if ((i <= nArg) && isempty(varargin{1}))
                nTR = size(varargin{1},1);
            end
            if nTR ~= length(TS{k}), error('Not all time series are the same length'); end
        end
        
        if (i<nArg) && iscell(varargin{i+1})
            tmp = varargin{i+1};
            if (length(tmp)==2 && iscell(tmp{2}))
                Titl{k} = tmp{1};
                Leg{k}  = tmp{2};
            else
                Leg{k} = tmp;
            end
            i = i+1;
        end
    end
    i = i+1;
end
nTS = k;
TS(k+1:nArg) = [];
Leg(k+1:nArg) = [];
Titl(k+1:nArg) = [];

return;


function h = MakeAxesPos(nTS,FigYsp)
%---------------------------------------------------------------------------
% MakeAxesPos: reserve space, and create the handle for time series plots.
% FORMAT h = MakeAxesPos(nTS, FigYsp)
% nTS        - number of time series.
% FigYsp     - reserved space
% h (output) - handles of the time series plots
%---------------------------------------------------------------------------
if isempty(FigYsp)
    FigYsp = [0 1];
elseif length(FigYsp) == 1
    FigYsp = [FigYsp 1];
end

VSpTop  = 0.05 + FigYsp(1);
VSpSpce = 0.03;
VSpBot  = 1+0.05-FigYsp(2);

VSpAxis = (1 - VSpTop - VSpBot - VSpSpce*(nTS-1))/nTS;

HSpLft  = 0.08;
HSpRgh  = 0.05;

HSpAxis = 1 - HSpLft - HSpRgh;
h = cell(nTS,1);

for i=1:nTS
    h{i} = [HSpLft  VSpBot+(nTS-i)*(VSpAxis+VSpSpce) HSpAxis VSpAxis];
end

spmd_pointer('Init');
return


function InitTimeCurs(nTR)
%---------------------------------------------------------------------------
% InitTimeCurs: set the initial temporal position.
% FORMAT InitTimeCurs(nTR)
% nTR    - number of scans
%---------------------------------------------------------------------------
global TimeCurs

if isempty(TimeCurs)
    TimeCurs = fix(nTR/2);
end
return


function SetTimeCurs(t)
%---------------------------------------------------------------------------
% SetTimeCurs: create the time cursor in the time series plots.
% FORMAT SetTimeCurs(t)
% t      - temporal position of the cursor to be set.
%---------------------------------------------------------------------------
global TimeCurs

if isempty(t)
    set(findobj('Tag','TimeCursor'),'visible','off');
    return;
else
    set(findobj('Tag','TimeCursor'),'Xdata',[t t]);
    set(findobj('Tag','TimeCursor'),'visible','on');
end

TimeCurs = t;
h = findobj('Tag','SPMd_MD');
s = get(h,'Name');
d = findstr(s,'Scan');


if ~isempty(d)
    s(d:end) = [];
else
    s = [s ' '];
end

s = [s sprintf('Scan %d',t)];
set(h,'Name',s);
h = findobj('Tag','SPMd_SS');
s = get(h,'Name');
d = findstr(s,'Scan');

if ~isempty(d)
    s(d:end) = [];
else
    s = [s ' '];
end

s = [s sprintf('Scan Summary: Scan %d',t)];
set(h,'Name',s);

if (~isempty(TimeCurs))
    spmd_pointer('Set',TimeCurs);
    return
end
return

function ClickMove
%---------------------------------------------------------------------------
% ClickMove: Move the cursor in other temporal viewers to the temporal
% position in the current time series plot.
%---------------------------------------------------------------------------
c   = get(gcbo,'CurrentPoint');
c   = round(c(1));
if (c<1), c = 1; end
spmd_mtsview('set',c)
return


function t = GetTimeCurs
%---------------------------------------------------------------------------
% GetTimeCurs: get the temporal position of the cursor.
%---------------------------------------------------------------------------
global TimeCurs
t = TimeCurs;
return


function UpdateTS(AxHnd,varargin)
%---------------------------------------------------------------------------
% Update: Update plot data w/out recreating axes, etc
% FORMAT UpdateTS(AxHnd,varargin)
%  AxHd       - Handles of axes to update
%  TS1,TS2... - Time series (Note, there must one for each row of AxHd)
%---------------------------------------------------------------------------
FS = spm('FontSizes');
[TSs,Titl,Legs,nTS,nTR] = countTS(varargin{:});
Tcurs = GetTimeCurs;

%- if TimeCurs does not exist, do nothing
if isempty(Tcurs)
    return;
end

XlimTR = [0.5 nTR+0.5];

for i = 1:nTS
    TS = TSs{i};
    AxHndChHnd  = [];
    
    Tmin = 100000;
    Tmax = 0;
    
    [numrow, numcol] = size(TS);
    
    %- Find maximum and minimum of TS values
    %-------------------------------------------------------------------------
    for j=1:numcol
        if Tmin > min(TS(:,j)) Tmin = min(TS(:,j)); end
        if Tmax < max(TS(:,j)) Tmax = max(TS(:,j)); end
    end
    
    if Tmin > 0
        if (Tmax - Tmin)> 50
            YlimTR = [0.99*Tmin 1.01*Tmax];        % for large scale
        elseif (Tmax - Tmin)>1
            YlimTR = [0.9*Tmin 1.1*Tmax];          % for small scale
        else
            YlimTR = [0.5*Tmin 1.5*Tmax];          % for tiny scale
        end
    else
        Yscale = max(abs(Tmin), Tmax);
        if abs(Tmin)<1
            YlimTR = [-(Yscale+0.2) (Yscale+0.2)];
        else
            YlimTR = [-Yscale*1.1 Yscale*1.1];
        end
    end
    
    
    axes(AxHnd(i));                           % specify current axes
    axis([XlimTR YlimTR]);                    % specify new axis in axes
    
    AxHndCh = get(AxHnd(i),'children');       % Get handle to the corresponding axes
    for j =1:length(AxHndCh)
        if strcmp(get(AxHndCh(j),'Tag'),'Temporal')
            AxHndChHnd = [AxHndChHnd,AxHndCh(j)]; % Get handles to the
            % 'temporal' plots
        end
    end
    
    AxHndChHnd = sort(AxHndChHnd);
    for j = 1:length(AxHndChHnd)
        set(AxHndChHnd(j),'Ydata',TS(:,j));     %-Update the TS figure
    end
    
    %- redraw cursor line
    %------------------------------------------------------------------------
    hCur = findobj(AxHnd(i),'Tag','TimeCursor');
    set(hCur, 'Xdata',Tcurs*[1 1],'Ydata',get(gca,'Ylim'),...
        'color','red','LineStyle',':',...
        'Tag','TimeCursor','HitTest','off');
    
    set(gca,'xlim',XlimTR,'UserData', hCur, 'fontsize', FS(8));
    
    %- redraw legend
    %------------------------------------------------------------------------
    if ~isempty(Legs{i})
        hl = legend(AxHndChHnd,Legs{i});
        set(hl,'fontsize',FS(8),'Tag','legend','HitTest','off');
    end
end
