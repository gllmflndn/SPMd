function varargout = spm_diagnosis(varargin)
% SPMd: Statistical Parametric Mapping Diagnosis (startup function)
% _________________________________________________________________________
% FORMAT spm_diagnosis('Ver')
% Show the name of toolbox
%
% FORMAT spm_diagnosis('Setup')
% Set up the SPMd main frame
%
% FORMAT spm_diagnosis('compui')
% Set up the menu for computation
%
% FORMAT spm_diagnosis('visualui')
% Set up the menu for visualization
%
% FORMAT spm_diagnosis('ms')
% Set up GUI for the model summary window: for both computation and
% visualization
%
% FORMAT spm_diagnosis('ss')
% Set up GUI for the scan summary window: for both computation and
% visualization
%
% FORMAT spm_diagnosis('sd')
% Set up GUI for the scan summary window: for both computation and
% visualization
%
% FORMAT spm_diagnosis('setimg')
% Draw results section buttons in Interactive window
%
% FORMAT spm_diagnosis('go')
% Manipulate the selection
%
% FORMAT spm_diagnosis('delete')
% Delete the objects with input tag or handles within
% the main menu figure
% _________________________________________________________________________
%
% SPMd is a toolbox for the diagnosis for the analysis functional brain
% mapping experiments using SPM.
%
% The spm_diagnosis function must startup in the directory containing
% SPM.mat file. Users need to run the SPM analysis for the data before
% using SPMd.
%
% This spm_diagnosis function displays a window with buttons leading to
% different aspects about diagnostic process. There are three buttons in
% the main frame:
% Compute: Compute all necessary spatial and temporal statistics that
%          are useful for the diagnosis
%          Three sub menu in the comput option:
%          Model summary - spatial
%          Scan summary  - temporal
%          Scan detail   - explore spatial detail
% Visualization: Display the choices of images of diagnosis statistics
%          Four sub menu in the visualization option:
%          Model summary - Display the choices of spatial images of
%                          diagnosis satatistics
%          Scan summary  - Display the choices of time series of
%                          diagnosis statistics
%          Model detail  - Display mean image and consecutive residual
%                          images around a specific time point
%          Scan detail   - Display various diagnosis plots at a specific
%                          voxel
%
% This diagnostic toolbox is designed for the linear regression model
% used in SPM. right now, it allows for models:
%          o Single session
%          o None, one or more conditions
%          o None or some additional covariates
%          o with/without high pass filtering
%          o with/without global scaling
% Only one serious restrictions for the use of SPMd:
%          o multiple sessions
% _________________________________________________________________________
%
% Reference:
% Luo, W-L and Nichols T. E. (2002) Diagnosis and Exploration of
% Massively Univariate fMRI Models. NeuroImage,19:1014-1032, 2003
%__________________________________________________________________________


if ~nargin, Action = 'SetUp'; else, Action = varargin{1}; end

%==========================================================================
switch lower(Action)
    
    %======================================================================
    case 'ver'                                            %-Version control
    %======================================================================
    
        varargout = {'SPMd'};
    
    %======================================================================
    case 'setup'                                     %-Set up the main menu
    %======================================================================
    
        %-Figure parameters
        %------------------------------------------------------------------
        FS = spm('FontSizes');
        PF = spm_platform('fonts');
        
        %-Create figure for SPMd menu
        %------------------------------------------------------------------
        RectW = [0 0 300 314];
        Rect0 = spm('WinSize','0',1);
        if all(ismember(Rect0(:),[0 1])), return; end
        Pos   = [Rect0(1)+(Rect0(3)-RectW(3))/2,...
            Rect0(2)+(Rect0(4)-RectW(4))/2,...
            RectW(3),...
            RectW(4)];
        Finter = figure('Color',[1 1 1]*.8,...
            'Name',spm_diagnosis('Ver'),...
            'NumberTitle','off',...
            'Position',Pos,...
            'Resize','on',...
            'Tag','Diag Menu',...
            'MenuBar','none',...
            'Visible','off');
        
        %-Set up buttons
        %------------------------------------------------------------------
        hReg = uicontrol(Finter,...
            'Style','Frame',...
            'Unit','normalized',...
            'Position',[0 0.9 1 0.1],...
            'BackgroundColor',spm('Colour'),...
            'Tag','Main');
        
        hReg = [hReg uicontrol(Finter,...
            'Style','PopUp',...
            'String','Compute|Model Summary|Scan Summary|Scan Detail',...
            'Unit','normalized',...
            'Position',[0.01 0.91 0.38 0.08],...
            'Callback','spm_diagnosis(''compUI'')',...
            'FontName',PF.times,...
            'FontWeight','Normal',...
            'FontAngle','Normal',...
            'FontSize',FS(12),...
            'HorizontalAlignment','center',...
            'Interruptible','on',...
            'Enable','on',...
            'ForegroundColor','b')];
        
        hReg = [hReg uicontrol(Finter,...
            'Style','PopUp',...
            'String','Visualization|Scan Summary|Model Summary|Model Detail|Scan Detail',...
            'Unit','normalized',...
            'Position',[0.4 0.91 0.38 0.08],...
            'Callback','spm_diagnosis(''visualUI'')',...
            'FontName',PF.times,...
            'FontWeight','Normal',...
            'FontAngle','normal',...
            'FontSize',FS(12),...
            'HorizontalAlignment','center',...
            'ForegroundColor','b')];
        
        hReg = [hReg uicontrol(Finter,...
            'Style','PushButton',...
            'String','exit',...
            'FontSize',FS(12),...
            'ForegroundColor','r',...
            'Callback','spm_diagnosis(''quit'')',...
            'Interruptible','on',...
            'Enable','on',...
            'Unit','normalized',...
            'Position',[0.79 0.91 0.2 0.08])];
        
        set(Finter,'Visible','on')
        set(hReg,'UserData',hReg);
        
        if nargout, varargout = {Finter}; end
        
    %======================================================================
    case 'compui'                                        %-Computation Menu
    %======================================================================
        %spm_diagnosis('compUI')
        
        %-Clear all previous sub-user interfaces
        %------------------------------------------------------------------
        Finter = spm_figure('findwin','Diag Menu');
        hMain  = get(findobj(Finter,'Tag','Main'),'userdata');
        hAll   = get(Finter,'Children');
        
        hdiff = setdiff(hAll,hMain');
        if ~isempty(hdiff)
            spm_diagnosis('delete','handle',hdiff)
        end
        
        %-Start to create the GUI menu
        %------------------------------------------------------------------
        hPM = gcbo;
        v = get(hPM,'Value');
        
        switch v
            case 1
                return
            case 2
                spm_diagnosis('ms',[],'compute');
            case 3
                spm_diagnosis('ss',[],'compute');
            case 4
                spm_diagnosis('sd',[],'compute');
            otherwise
                error('Unknown computation option.');
        end
        
        set(hPM,'Value',1);
        
    %======================================================================
    case 'visualui'                                    %-Visualization Menu
    %======================================================================
        %spm_diagnosis('visualUI')

        %-Clear all previous sub-user interfaces
        %------------------------------------------------------------------
        Finter = spm_figure('findwin','Diag Menu');
        hMain  = get(findobj(Finter,'Tag','Main'),'userdata');
        hAll   = get(Finter,'Children');
        
        hdiff = setdiff(hAll,hMain');
        spm_diagnosis('delete','handle',hdiff)
        
        %-Start to create the GUI menu
        %------------------------------------------------------------------
        hPM = gcbo;
        v = get(hPM,'Value');
        
        switch v
            case 1
                return
            case 2
                spm_diagnosis('ss',[],'visual');
            case 3
                spm_diagnosis('ms',[],'visual');
            case 4
                spmd_MD;
            case 5
                spmd_SD;
            otherwise
                error('Unknown visualization option.');
        end
        
        set(hPM,'Value',1);
        
    %======================================================================
    case 'ms'          %-Model summary menu (computation and visualization)
    %======================================================================
        %spm_diagnosis('ms',Finter,type)
        
        %-Figure parameters
        %------------------------------------------------------------------
        WS     = spm('WinScale');
        FS     = spm('FontSizes');
        PF     = spm_platform('fonts');
        
        if nargin < 3
            type = 'compute';
        else
            type = varargin{3};
        end
        
        if nargin < 2 || isempty(varargin{2})
            Finter = 'Diag Menu';
        else
            Finter = varargin{2};
        end
        Finter = spm_figure('Getwin',Finter);
        
        %-Load SPM.mat
        %------------------------------------------------------------------
        if exist(fullfile(pwd,'SPMd_MS.mat'),'file')==2
            load SPMd_MS MS
            swd = MS.SPM.swd;
        elseif exist(fullfile(pwd,'SPMd_SS.mat'),'file')==2
            load SPMd_SS SS
            swd = SS.SPM.swd;
        else
            swd = pwd;
        end
        
        if exist(fullfile(swd,'SPM.mat'),'file')~=2
            [spmmat_handle, sts] = spm_select(1,'^SPM\.mat$','Check whether you have SPM.mat file for analysis');
            if ~sts
                error('No SPM.mat in current directory.');
            end
            swd = spm_file(spmmat_handle,'fpath');
        end
        
        xSPM = load(fullfile(swd,'SPM.mat'),'SPM');
        xSPM = xSPM.SPM;
        
        dimX = size(xSPM.xX.X,2);
        
        %-Require parameters for the GUI figure
        %------------------------------------------------------------------
        D	  = get_MS_imgs;
        if strcmp(type,'compute')
            D = D(5:end);
            s = 1;
        elseif strcmp(type,'visual')
            D = D;
            s = 2;
        end
        
        nD	  = length(D);
        mD	  = ceil((nD+1)/3);
        Lx    = 0.01;
        LX    = 0.32;
        Ly    = 0.01;
        LY    = 0.05;
        
        %-Create frame for Results GUI objects
        %------------------------------------------------------------------
        Ypos      = 0.84-((Ly+LY)*mD+Ly);
        Ylength   = (LY+Ly)*mD+Ly;
        hMS_UI   = uicontrol(Finter,...
            'Style','Frame',...
            'Tag','hMS_UI',...
            'Unit','normalized',...
            'Position',[0.00 Ypos 1.00 Ylength],...
            'BackgroundColor',spm('Colour'));
        
        %-Create the text at the top
        %------------------------------------------------------------------
        h = [];
        h = [h uicontrol(Finter,...
            'Style','Text',...
            'String','Select Model Summary Images',...
            'Unit','normalized',...
            'Position',[0.0 0.84 1 0.06],...
            'FontName',PF.times,...
            'FontWeight','Bold',...
            'FontAngle','Italic',...
            'FontSize',FS(14),...
            'HorizontalAlignment','left',...
            'ForegroundColor','b')];
        
        %-Buttons to select spatial statistical images
        %------------------------------------------------------------------
        Ch	= [];
        
        if strcmp(type,'visual')
            Ch	= [Ch uicontrol(Finter,...
                'Style','checkbox',...
                'String',D(1).label,...
                'Min', 0, 'Max',4,...
                'Value',0,...
                'CallBack','spm_diagnosis(''SetImg'',''test'')',...
                'Unit','normalized',...
                'Position',[Lx 0.84-1*(Ly+LY) LX LY],...
                'Interruptible','on',...
                'ToolTipString',D(1).help,...
                'UserData',D(1),...
                'ForegroundColor','k')];
            
            Ch	= [Ch uicontrol(Finter,...
                'Style','checkbox',...
                'String',D(2).label,...
                'Min', 0, 'Max',4,...
                'Value',0,...
                'CallBack','spm_diagnosis(''SetImg'',''con'')',...
                'Unit','normalized',...
                'Position',[Lx+(LX+Lx) 0.84-1*(Ly+LY) LX LY],...
                'Interruptible','on',...
                'UserData',D(2),...
                'ForegroundColor','k')];
            
            Ch	= [Ch uicontrol(Finter,...
                'Style','checkbox',...
                'String',D(3).label,...
                'Min', 0, 'Max',4,...
                'Value',4,...
                'Unit','normalized',...
                'Position',[Lx+2*(LX+Lx) 0.84-1*(Ly+LY) LX LY],...
                'Interruptible','on',...
                'UserData',D(3),...
                'ForegroundColor','k')];
        end
        
        for i = 1:3
            for j = s:mD
                if 3*(j-1)+i < nD+1
                    
                    if (strcmp(D(3*(j-1)+i).label,'Homo. vs h(Y)') || ...
                            strcmp(D(3*(j-1)+i).label,'Homo. vs X'))
                        if (dimX==1 && strcmp(type,'visual'))
                            Ch = [Ch uicontrol(Finter,...
                                'Style','checkbox',...
                                'String',D(3*(j-1)+i).label,...
                                'Min', 0, 'Max',4,...
                                'Value',0,...
                                'HorizontalAlignment','Center',...
                                'Unit','normalized',...
                                'Position',[Lx+(i-1)*(LX+Lx) 0.84-(LY+Ly)*j LX LY],...
                                'ToolTipString',D(3*(j-1)+i).help,...
                                'UserData',D(3*(j-1)+i),...
                                'ForegroundColor','k')];
                        else
                            Ch = [Ch uicontrol(Finter,...
                                'Style','checkbox',...
                                'String',D(3*(j-1)+i).label,...
                                'Min', 0, 'Max',4,...
                                'Value',4,...
                                'HorizontalAlignment','Center',...
                                'Unit','normalized',...
                                'Position',[Lx+(i-1)*(LX+Lx) 0.84-(LY+Ly)*j LX LY],...
                                'ToolTipString',D(3*(j-1)+i).help,...
                                'UserData',D(3*(j-1)+i),...
                                'ForegroundColor','k')];
                        end
                    else
                        Ch = [Ch uicontrol(Finter,...
                            'Style','checkbox',...
                            'String',D(3*(j-1)+i).label,...
                            'Min', 0, 'Max',4,...
                            'Value',4,...
                            'HorizontalAlignment','Center',...
                            'Unit','normalized',...
                            'Position',[Lx+(i-1)*(LX+Lx) 0.84-(LY+Ly)*j LX LY],...
                            'ToolTipString',D(3*(j-1)+i).help,...
                            'UserData',D(3*(j-1)+i),...
                            'ForegroundColor','k')];
                    end
                end
            end
        end
        
        h = [h uicontrol(Finter,...
            'String','Cancel',...
            'Unit','normalized',...
            'Position',[Lx+2*(LX+Lx)+0.16 Ypos+Ly 0.16 LY],...
            'CallBack','spm_diagnosis(''delete'',''tag'',''hMS_UI'')',...
            'ForegroundColor','r')];
        
        %save selected choices
        if strcmp(type,'compute')
            h = [h uicontrol(Finter,...
                'String','Run',...
                'Unit','normalized',...
                'Position',[Lx+2*(Lx+LX) Ypos+Ly 0.16 LY],...
                'CallBack','spm_diagnosis(''go'',''ms_comp'');spm_diagnosis(''delete'',''tag'',''hMS_UI'')',...
                'ForegroundColor','r')];
            MS.Comp = Ch;
        elseif strcmp(type,'visual')
            h = [h uicontrol(Finter,...
                'String','Run',...
                'Unit','normalized',...
                'Position',[Lx+2*(Lx+LX) Ypos+Ly 0.16  LY],...
                'CallBack','spm_diagnosis(''go'',''ms_visual'');spm_diagnosis(''delete'',''tag'',''hMS_UI'')',...
                'ForegroundColor','r')];
            
            MS.Visual = Ch;
        end
        
        set(Finter,'UserData',MS)
        set(hMS_UI,'userdata',[hMS_UI h Ch])
        
    %======================================================================
    case 'ss'           %-Scan summary menu (computation and visualization)
    %======================================================================
        %spm_diagnosis('ss',Finter,type)
        
        %-window Parameters
        %------------------------------------------------------------------
        WS     = spm('WinScale');
        FS     = spm('FontSizes');
        PF     = spm_platform('fonts');
        
        if nargin < 3
            type = 'compute';
        else
            type = varargin{3};
        end
        
        if nargin < 2 || isempty(varargin{2})
            Finter = 'Diag Menu';
        else
            Finter = varargin{2};
        end
        Finter = spm_figure('Getwin',Finter);
        
        
        %-Require parameters for the GUI figure
        %------------------------------------------------------------------
        D       = get_SS_ts(type);
        nD	    = length(D);
        mD	    = ceil((nD+1)/2);
        Lx      = 0.02;
        LX      = 0.47;
        Ly      = 0.01;
        LY      = 0.06;
        
        %-Create frame for Results GUI objects
        %------------------------------------------------------------------
        Ypos      = 0.84-((Ly+LY)*mD+Ly);
        Ylength   = (LY+Ly)*mD+Ly;
        hSS_UI    = uicontrol(Finter,'Style','Frame','Tag','hSS_UI',...
            'unit','normalized',...
            'Position',[0.00 Ypos 1 Ylength],...
            'BackgroundColor',spm('Colour'));
        
        %-Create the text at the top
        %------------------------------------------------------------------
        %delete(h)
        h = [];
        h = [h uicontrol(Finter,...
            'Style','Text',...
            'String','Select Scan Summary Series',...
            'Unit','normalized',...
            'Position',[0 0.84 1 0.06],...
            'FontName',PF.times,...
            'FontWeight','Bold',...
            'FontAngle','Italic',...
            'FontSize',FS(14),...
            'HorizontalAlignment','left',...
            'ForegroundColor','b')];
        
        
        %-Buttons to select spatial statistical images
        %------------------------------------------------------------------
        %delete(Ch)
        Ch	= [];
        if (strcmp(type,'visual'))
            load SPMd_SS;
        end
        for i = 1:2
            for j = 1:mD
                if 2*(j-1)+i < nD+1
                    val=4;
                    
                    %-deselect Rotation and Shift option when unavailable
                    if (strcmp(D(2*(j-1)+i).label,'Rotation Parameters') ||...
                            strcmp(D(2*(j-1)+i).label,'Shift Parameters')) && ...
                            (isfield(SS,'RX') && isempty(SS.RX) && strcmp(type,'visual'))
                        val = 0;
                    end
                    
                    Ch = [Ch uicontrol(Finter,...
                        'Style','checkbox',...
                        'String',D(2*(j-1)+i).label,...
                        'Min', 0, 'Max',4,...
                        'value',val,...
                        'HorizontalAlignment','Center',...
                        'Unit','normalized',...
                        'Position',[Lx+(i-1)*(LX+Lx) 0.84-(LY+Ly)*j LX LY],...
                        'ToolTipString',D(2*(j-1)+i).help,...
                        'UserData',D(2*(j-1)+i),...
                        'ForegroundColor','k')];
                end
            end
        end
        
        h = [h uicontrol(Finter,...
            'String','Cancel',...
            'Unit','normalized',...
            'Position',[Lx+(LX+Lx)+0.25 Ypos+Ly 0.22 LY],...
            'CallBack','spm_diagnosis(''delete'',''tag'',''hSS_UI'')',...
            'ForegroundColor','r')];
        
        %save selected choices
        if strcmp(type,'compute')
            h = [h uicontrol(Finter,...
                'String','Run',...
                'Unit','normalized',...
                'Position',[Lx+(Lx+LX) Ypos+Ly 0.22 LY],...
                'CallBack','spm_diagnosis(''go'',''ss_comp'');spm_diagnosis(''delete'',''tag'',''hSS_UI'')',...
                'ForegroundColor','r')];
            SS.Comp = Ch;
        elseif strcmp(type,'visual')
            h = [h uicontrol(Finter,...
                'String','Run',...
                'Unit','normalized',...
                'Position',[Lx+(Lx+LX) Ypos+Ly 0.22 LY],...
                'CallBack','spm_diagnosis(''go'',''ss_visual'');spm_diagnosis(''delete'',''tag'',''hSS_UI'')',...
                'ForegroundColor','r')];
            SS.Visual = Ch;
        end
        
        set(Finter,'UserData',SS)
        set(hSS_UI,'userdata',[hSS_UI h Ch])
        
    %======================================================================
    case 'sd'            %-Scan detail menu (computation and visualization)
    %======================================================================
        %spm_diagnosis('sd',Finter,type)
        
        if nargin < 3
            type = 'compute';
        else
            type = varargin{3};
        end
        
        if nargin < 2 || isempty(varargin{2})
            Finter = 'Diag Menu';
        else
            Finter = varargin{2};
        end
        Finter = spm_figure('Getwin',Finter);
        
        if strcmp(type,'compute')
            
            %-Find the SPM.mat
            if exist(fullfile('.','SPMd_MS.mat'),'file')==2
                load SPMd_MS
                swd = MS.SPM.swd;
            elseif exist(fullfile('.','SPMd_SS.mat'),'file')==2
                load SPMd_SS
                swd = SS.SPM.swd;
            else
                swd = spm_str_manip(spm_select(1,'SPM.mat','Select SPM.mat'),'H');
            end
            
            if exist(fullfile(swd,'SPM.mat'),'file')~=2
                spmmat_handle = spm_select(1,'SPM.mat','Select SPM.mat file');
                swdSD = spm_str_manip(spmmat_handle,'hd');
                xSPM = load(fullfile(swdSD,'SPM.mat'));
            else
                xSPM = load(fullfile(swd,'SPM.mat'));
            end
            
            xSPM = xSPM.SPM;
            fprintf('\nPreparing spatiotemporal statistics... \n')
            PstRes = spmd_MkResid(xSPM);
            save SPMd_SD PstRes;
            fprintf('\nDone!')
            
        elseif strcmp(type,'visual')
            spmd_SD
        end
        
    %======================================================================
    case 'setimg'      %-Draw results section buttons in Interactive window
    %======================================================================
        % spm_diagnosis('SetImg',varagin)
        
        Ch	= get(findobj('Tag','Diag Menu'),'UserData');
        BaseNm	= varargin{1};
        
        switch varargin{2}
            
            case 'test'
                D	= get(Ch.Visual(1),'UserData');
                if get(Ch.Visual(1),'value')
                    D.fname	= spm_select(1,'spm.*\.img');
                else
                    D.fname	= '';
                end
                set(Ch.Visual(1),'UserData',D);
                
            case 'con'
                D	= get(Ch.Visual(2),'UserData');
                if get(Ch.Visual(2),'value')
                    D.fname	= spm_select(1,'con.*\.img');
                else
                    D.fname	= '';
                end
                set(Ch.Visual(2),'UserData',D);
                
        end
        
    %======================================================================
    case 'go'                                    %-Manipulate the selection
    %======================================================================
        %spm_diagnosis('go',type)
        
        type = varargin{2};
        
        switch type
            case 'ms_comp'
                Ch	= get(findobj('Tag','Diag Menu'),'UserData');
                Ch    = Ch.Comp;
                D	    = get_MS_imgs;
                nD	= length(D);
                Im	= {};
                
                for i = 1:nD-4
                    if get(Ch(i),'value')
                        DI	= get(Ch(i),'UserData');
                        Im	= {Im{:} DI.name};
                    end
                end
                
                if ~isempty(Im)
                    swd     = spm_str_manip(spm_select(1,'SPM.mat','Select SPM.mat'),'H');
                    xSPM    = load(fullfile(swd,'SPM.mat'));
                    if (length(fieldnames(xSPM)) == 1)
                        xSPM= xSPM.SPM;
                    end
                    spmd_comp_MS(xSPM,Im);
                end
                
                
            case 'ms_visual'
                Ch	= get(findobj('Tag','Diag Menu'),'UserData');
                Ch    = Ch.Visual;
                D	    = get_MS_imgs;
                nD	= length(D);
                Im	= {};
                Ds	= {};
                
                for i = 1:nD
                    if get(Ch(i),'value')
                        DI	= get(Ch(i),'UserData');
                        Im	= {Im{:},DI.fname};
                        Ds	= {Ds{:},DI.desc};
                    end
                end
                
                if isempty(Im)
                    spmd_MS;
                end
                
                spmd_MS('Img',1:length(Im),Im,'Desc',1:length(Ds),Ds);
                
            case 'ss_comp'
                Ch	= get(findobj('Tag','Diag Menu'),'UserData');
                Ch    = Ch.Comp;
                D	= get_SS_ts('compute');
                nD	= length(D);
                Fig	= {};
                
                for i = 1:nD
                    if get(Ch(i),'value')
                        DI	= get(Ch(i),'UserData');
                        Fig	= {Fig{:},DI.fname};
                    end
                end
                
                nFig	= length(Fig);
                spmd_comp_SS(Fig{1:nFig});
                
            case 'ss_visual'
                Ch	= get(findobj('Tag','Diag Menu'),'UserData');
                Ch    = Ch.Visual;
                D	= get_SS_ts('visual');
                nD	= length(D);
                Fig	= {};
                
                for i = 1:nD
                    if get(Ch(i),'value')
                        DI	= get(Ch(i),'UserData');
                        Fig	= {Fig{:},DI.fname};
                    end
                end
                
                if isempty(Fig)
                    spmd_SS
                end
                
                nFig	= length(Fig);
                spmd_SS('ts',{Fig{1:nFig}});
                
            otherwise
                error('Unknown command');
        end
        
    %======================================================================
    case 'delete'              %-Delete objects within the main menu figure
    %======================================================================
        %spm_diagnosis('delete',type,var)
        %  type = 'handle'; var = handle;
        %  type = 'tag';    var = tag;
        
        if nargin < 2, return, end
        
        type = varargin{2};
        
        switch type
            case 'handle'
                C   = varargin{3};
            case 'tag'
                tag = varargin{3};
                F   = spm_figure('FindWin','Diag Menu');
                h   = findobj(F,'Tag',tag);
                C   = get(h,'UserData');
        end
        
        spm_results_ui('Delete',C);
        
    %======================================================================
    case 'quit'                                                %-Leave SPMd
    %======================================================================
        %spm_diagnosis('Quit')
    
        Finter = spm_figure('FindWin','Diag Menu');
        close(Finter);
        
    %======================================================================
    otherwise
    %======================================================================
        error(['Unknown SPMd action (' Action ')']);
        
end


function D = get_MS_imgs
%==========================================================================
% Build the menu structure of all diagnostic and exploratory statistics
%==========================================================================

D(1).fname	= '';
D(1).desc	= 'Test statistic image';
D(1).name 	= 'Test Statistic';
D(1).label  = 'Test Stat';
D(1).help	= 'Hypothesis testing statistics of interest';

D(2).fname	= '';
D(2).desc	= 'Contrast image';
D(2).name	= 'Contrast';
D(2).label  = 'Contrast';
D(2).help	= 'Contrast Image of interest';

D(3).fname	= 'Mean.nii';
D(3).desc	= 'Mean image';
D(3).name	= 'Mean';
D(3).label  = 'Mean';
D(3).help	= 'Mean Image';

D(4).fname	= 'SPMd_ResRMS.nii';
D(4).desc	= 'Stdev of residuals';
D(4).name	= 'ResRMS';
D(4).label  = 'ResRMS';
D(4).help	= 'Standard deviation of residuals';

D(5).fname	= 'SPMd_PCorr.nii';
D(5).desc	= '-log_{10}p DW';
D(5).name	= 'Corr';
D(5).label  = 'Correlation';
D(5).help	= 'Durbin-Watson statistic assessing the autocorrelation of residuals';

D(6).fname	= 'SPMd_PDep.nii';
D(6).desc	= '-log_{10}p CP';
D(6).name	= 'Dep';
D(6).label  = 'Dependence';
D(6).help	= 'Cumulative periodogram assessing white noise assumption of residuals';

D(7).fname	= 'SPMd_PHomo1.nii';
D(7).desc	= '-log_{10}p CW: Global';
D(7).name	= 'Homo1';
D(7).label  = 'Homo. vs Glob.';
D(7).help	= ['Cook-Weisberg score statistic assessing constant variance assumption ' ...
    'of residuals wrt residuals wrt global signal'];

D(8).fname	= 'SPMd_PHomo2.nii';
D(8).desc	= '-log_{10}p CW: Predicted';
D(8).name	= 'Homo2';
D(8).label  = 'Homo. vs h(Y)';
D(8).help	= ['Cook-Weisberg score statistic assessing constant variance assumption ' ...
    'of residuals wrt predicted response'];

D(9).fname	= 'SPMd_PHomo3.nii';
D(9).desc	= '-log_{10}p CW: Active';
D(9).name 	= 'Homo3';
D(9).label  = 'Homo. vs X';
D(9).help	= ['Cook-Weisberg score statistic assessing constant variance assumption ' ...
    'of residuals wrt experimental condition'];

D(10).fname	= 'SPMd_PNorm.nii';
D(10).desc	= '-log_{10}p SW';
D(10).name	= 'Norm';
D(10).label = 'Normality';
D(10).help	= 'Shapiro-Wilk statistic assessing normality assumption of residuals';

D(11).fname	= 'SPMd_Outl.nii';
D(11).desc	= 'Outlier count';
D(11).name	= 'Outl';
D(11).label = 'Outlier';
D(11).help	= 'Assess the number of outliers at the voxel';


function D = get_SS_ts(type)
%==========================================================================
% Build the menu structure of all diagnostic and exploratory statistics
%==========================================================================
D(1).fname	= 'predint';
D(1).label      = 'Predictor';
D(1).help	= 'The predicted signal';

D(2).fname	= 'global';
D(2).label      = 'Global Signal';
D(2).help	= 'The temporal global intensity';

D(3).fname	= 'toutlier';
D(3).label      = 'Outlier Rate';
D(3).help	= 'Percent of expected outliers';

if strcmp(type,'compute')
    D(4).fname	= 'regparm';
    D(4).label    = 'Reg. Parameters';
    D(4).help	= 'Registration shift and rotation movement parameters';
    
    D(5).fname	= 'pg';
    D(5).label    = 'Average periodogram';
    D(5).help	= ['Average periodogram of raw residuals of',...
        'the time series over the brain'];
elseif strcmp(type,'visual')
    D(4).fname	= 'shift';
    D(4).label    = 'Shift Parameters';
    D(4).help	= 'Registration shift movement parameters';
    
    D(5).fname	= 'rotate';
    D(5).label    = 'Rotation Parameters';
    D(5).help	= 'Registration rotation movement parameters';
    
    D(6).fname	= 'pg';
    D(6).label    = 'Average periodogram';
    D(6).help	= ['Average periodogram of raw residuals of',...
        'the time series over the brain'];
end
