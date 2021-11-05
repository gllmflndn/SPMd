function spmd_SptlBar(varargin)
% Create toolbar in the figure
% FORMAT spmd_SptlBar(varargin)
%   varargin - different types of options
% _________________________________________________________________________
% @(#)spmd_SptlBar.m    1.1 04/06/22

% _______________________ Functions called ________________________________
%
%       spm_figure
%       spmd_orthviews
% _________________________________________________________________________


%-condition arguments
%--------------------------------------------------------------------------
if (nargin == 0)
    Action = 'Create';
else
    Action = varargin{1};
end


switch lower(Action)
    case 'create'
        %==================================================================
        % spmd_SptlBar('create',fig)
        
        if nargin<2
            if ~isempty(get(0,'Children'))
                F=gcf;
            else
                F='';
            end
            
        else
            F=fig;
        end
        
        F = spm_figure('FindWin',F);
        if isempty(F)
            return,
        end
        
        %-Get position and size parameters
        %------------------------------------------------------------------
        cUnits = get(F,'Units');
        set(F,'Units','Pixels');
        P     = get(F,'Position'); P  = P(3:4);     % Figure dimensions {pixels}
        S_Gra = P./[600, 865];              % x & y scaling coefs
        
        nBut  = 6;
        nGap  = 2;
        sx    = floor(P(1)./(nBut+(nGap+2)/6));     % uicontrol object width
        dx    = floor(2*sx/6);              % inter-uicontrol gap
        sy    = floor(20*S_Gra(1));         % uicontrol object height
        x0    = dx;                 % initial x position
        x     = dx/2;                   % uicontrol x position
        y     = P(2) - sy;              % uicontrol y position
        y2    = P(2) - 2.25*sy;             % uicontrol y position
        FS    = round(10*min(S_Gra));           % uicontrol font size
        
        %-Delete any existing 'ToolBar' 'Tag'ged objects
        %------------------------------------------------------------------
        cSHH = get(0,'ShowHiddenHandles');
        set(0,'ShowHiddenHandles','on')
        delete(findobj(F,'Tag','ToolBar'));
        set(0,'ShowHiddenHandles',cSHH)
        
        %-Create Frame for controls
        %------------------------------------------------------------------
        uipanel(F,...
            'Position',[-4 (P(2) - 1.25*sy) P(1)+8 1.25*sy+4],...
            'Tag','ToolBar',...
            'HandleVisibility','callback');
        
        %-Create uicontrol objects
        %------------------------------------------------------------------
        uicontrol(F,'String','Print','ToolTipString','print figure',...
            'Position',[x y sx sy],...
            'CallBack','spm_figure(''Print'',gcf)',...
            'FontSize',FS,...
            'Interruptible','off','BusyAction','cancel',...
            'Tag','ToolBar','HandleVisibility','callback',...
            'ForegroundColor','b'); x = x+sx;
        
        uicontrol(F,'String','Move','ToolTipString',...
            'Move the crosshair in corresponding figure windows',...
            'Position',[x y sx sy],...
            'CallBack','spmd_SptlBar(''FindFig'')',...
            'FontSize',FS,...
            'Interruptible','off','BusyAction','cancel',...
            'Tag','ToolBar','HandleVisibility','callback',...
            'ForegroundColor','b'); x = x+sx;
        
        uicontrol(F,'String','Switch',...
            'ToolTipString','Switch the spatial image windows',...
            'Position',[x y sx sy],...
            'CallBack','spmd_SptlBar(''Switch'')',...
            'FontSize',FS,...
            'Interruptible','off','BusyAction','cancel',...
            'Tag','ToolBar','HandleVisibility','callback',...
            'ForegroundColor','b'); x = x+sx;
        
        set(F,'Units',cUnits)
        
        
    case 'findfig'
        %==================================================================
        % spmd_SptlBar('findfig')
        
        fig0 = findobj('Tag','Graphics');
        fig1 = findobj('Tag','Model Summary');
        fig2 = findobj('Tag','Scan Detail');
        
        figs = [fig0 fig1 fig2];
        
        spmd_orthviews('fig',gcf,'reposother',figs);
        
    case 'switch'
        %==================================================================
        % spmd_SptlBar('switch')
        h = gcf;
        
        fig0 = findobj('Tag','Graphics');
        fig1 = findobj('Tag','Model Summary');
        fig2 = findobj('Tag','Scan Detail');
        
        figs = [fig0 fig1 fig2];
        
        for i = 1:length(figs)
            set(figs(i),'position',get(h,'position'));
            if figs(i) == h
                set(figs(i),'visible','off');
            else
                set(figs(i),'visible','on');
            end
        end
        
end
