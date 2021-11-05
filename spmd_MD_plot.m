function spmd_MD_plot(h,type,D)
% Create and update residual plots for the current spatial voxel
% FORMAT spmd_MD_plot(h,type,D)
% h    - current axes' handle
% type - users selected plots
%        o Residual vs experimental predictor
%        o Residual vs global signal
%        o Absolute value of residual vs. experimental predictor
%        o Absolute value of residual vs. global signal
%        o Lag-1 serial plot
%        o Normality plot
% D    - all necessary data for the plot
%        o Residual
%        o Absolute value of residuals
%        o Studentized residuals
%        o Raw data
%        o Fitted drift data in the fMRI time series
%        o Fitted data of full model
%        o Global Signal
% _________________________________________________________________________
% @(#)spmd_MD_plot.m

% ______________________________Function Called __________________________
%
%        spm_platform
%        spmd_Nplot
%        plot_res (internal)
% ________________________________________________________________________

if nargin < 3
    if(iscell(h))
        h = celldisp(h);
        D = get(h,'userdata');
    else
        D = get(h,'userdata');
    end
    hf = 0;
else
    hf = 1;
end
if nargin < 2, type = 1; end

%cla(h)
axes(h);
switch lower(type)
    case {1,'exp_r'}
        X      = D(7);
        Y      = D(1);
        titl   = 'Res vs Experimental Predictor';
        plot_res(h,X,Y,titl,D, hf);
        
    case {2,'global_r'}
        X      = D(8);
        Y      = D(1);
        titl   = 'Res vs Global Signal';
        plot_res(h,X,Y,titl,D, hf);
        
    case {3,'exp_abs'}
        X      = D(7);
        Y      = D(2);
        titl   = '|Res| vs Experimental Predictor';
        plot_res(h,X,Y,titl,D, hf);
        
    case {4,'global_abs'}
        X      = D(8);
        Y      = D(2);
        titl   = '|Res| vs Global Signal';
        plot_res(h,X,Y,titl,D, hf);
        
    case {5,'lag_1'}
        Y.data = [D(1).data(1:end-1);0];
        Y.desp = 'Residual(i)';
        Y.help = 'Residual at time point i';
        Y.abline = 'h=0';
        X.data = [D(1).data(2:end);0];
        X.desp = 'Residual(i-1)';
        X.help = 'Residual at time point i-1';
        X.abline = 'h=0';
        titl   = 'Lag-1 plot';
        plot_res(h,X,Y,titl,D, hf);
        
    case {6,'QQ'}
        Y      = D(1).data;
        spmd_Nplot(Y,h,D);
end


function plot_res(h,X,Y,titl,D, hf)
%==========================================================================
% redraw the current plot
% h - current axes' handler
% X - xdata
% Y - ydata
% D - all data information set to this plot
%==========================================================================

global TimeCurs
%-window Parameters
%--------------------------------------------------------------------------
WS     = spm('WinScale');
FS     = spm('FontSizes');
PF     = spm_platform('fonts');

Xdata = X.data;
Ydata = Y.data;

if hf==0
    hp = plot(Xdata,Ydata,'.');   %- new plot
    set(hp,'Tag', 'MD_pointer');
    if strcmp(X.abline,'h=0')
        hy = get(h,'Ylim');
        %- draw a vertical line
        lv = line([0 0], hy, 'linestyle','-.','color','r','Tag', 'vt=0');
    end
    
    if strcmp(Y.abline,'h=0')
        hx = get(h,'Xlim');
        %- draw a horizontal line
        lz = line(hx, [0 0], 'linestyle','-.','color','r','Tag', 'hz=0');
    end
    
    hold on;
    %- plot a red star
    ht = plot(Xdata(TimeCurs), Ydata(TimeCurs),'r*', 'Tag', 'SelPt');
    hold off;
    Ux = uicontextmenu;
    hux = xlabel(X.desp,'uicontextmenu',Ux);
    set(hux,'fontsize',FS(8),'Tag','xaxlable');
    
    Uy = uicontextmenu;
    huy = ylabel(Y.desp,'uicontextmenu',Uy);
    set(huy,'fontsize',FS(8),'Tag', 'yaxlable');
    
    xdesp = uimenu(Ux,'Label',X.help);
    ydesp = uimenu(Uy,'Label',Y.help);
    
    htt = title(titl,'uicontextmenu',Ux);
    set(htt,'fontsize',FS(8),'Tag','ttitle');
    
    
else
    
    hp = findobj(h,'Tag','MD_pointer');   %- find MD_pointer
    if (size(Xdata,2) == 1)
        set(hp, 'Xdata', Xdata, 'Ydata', Ydata, 'Tag', 'MD_pointer'); %- update
        
        %- adjuste 'xlim' and 'ylim'
        
        if min(Xdata)> 1000
            minXdata = min(Xdata)*0.9995;
            maxXdata = max(Xdata)*1.0005;
        elseif (min(Xdata)<1000) & (min(Xdata)>=100)
            minXdata = min(Xdata)*0.995;
            maxXdata = max(Xdata)*1.005;
        elseif (min(Xdata)<100) & (min(Xdata) >=10)
            minXdata = min(Xdata)*0.95;
            maxXdata = max(Xdata)*1.05;
        elseif (min(Xdata)<10) & (min(Xdata) >= 1)
            minXdata = min(Xdata)*0.5;
            maxXdata = max(Xdata)*1.5;
        elseif (min(Xdata)<1) & (min(Xdata)>=0)
            minXdata = min(Xdata)*0.5;
            maxXdata = max(Xdata)*1.5;
        elseif (min(Xdata)<0) & (min(Xdata)>=-1)
            minXdata = min(Xdata)*2;
            maxXdata = max(Xdata)*2;
        elseif (min(Xdata)<-1) & (min(Xdata)>=-10)
            minXdata = min(Xdata)*1.5;
            maxXdata = max(Xdata)*1.5;
        else
            minXdata = min(Xdata)*1.1;
            maxXdata = max(Xdata)*1.1;
        end
        
        if min(Ydata) >1000
            minYdata = min(Ydata)*0.9995;
            maxYdata = max(Ydata)*1.0005;
        elseif (min(Ydata)<1000) & (min(Ydata)>=100)
            minYdata = min(Ydata)*0.995;
            maxYdata = max(Ydata)*1.005;
        elseif (min(Ydata)<100) & (min(Ydata)>=10)
            minYdata = min(Ydata)*0.95;
            maxYdata = max(Ydata)*1.05;
        elseif (min(Ydata)<10) & (min(Ydata)>=1)
            minYdata = min(Ydata)*0.95;
            maxYdata = max(Ydata)*1.05;
        elseif (min(Ydata)<1) & (min(Ydata)>=0)
            minYdata = min(Ydata)*0.5;
            maxYdata = max(Ydata)*1.5;
        elseif (min(Ydata)<0) & (min(Ydata)>=-1)
            minYdata = min(Ydata)*2;
            maxYdata = max(Ydata)*2;
        elseif (min(Ydata)<-1) & (min(Ydata)>-10)
            minYdata = min(Ydata)*1.5;
            maxYdata = max(Ydata)*1.5;
        else
            minYdata = min(Ydata)*1.1;
            maxYdata = max(Ydata)*1.1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Tom's change
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (minXdata == maxXdata)
            minXdata = minXdata - 0.5;
            maxXdata = maxXdata + 0.5;
        end
        if (minYdata == maxYdata)
            minYdata = minYdata - 0.5;
            maxYdata = maxYdata + 0.5;
        end
        
        
        %- set new axis
        set(h,'xlim',[minXdata maxXdata] ,'ylim', [minYdata, maxYdata]);
        
        lz = findobj(h, 'Tag', 'hz=0');
        lv = findobj(h, 'Tag', 'vt=0');
        ht = findobj(h, 'Tag', 'SelPt');
        Ux = get(h, 'uicontextmenu');
        
        if strcmp(X.abline,'h=0')
            hy = get(h,'Ylim');
            if ~isempty(lv)
                %- update vertical line information
                set(lv, 'Xdata', [0 0], 'Ydata', hy, 'linestyle','-.','color','r', 'Tag', 'vt=0');
            else
                lv = line([0 0], hy, 'linestyle','-.','color','r','Tag', 'vt=0');
            end
        end
        
        if strcmp(Y.abline,'h=0')
            hx = get(h,'Xlim');
            if ~isempty(lz)
                %- update vertical line information
                set(lz, 'Xdata', hx, 'Ydata', [0 0], 'linestyle','-.','color','r', 'Tag', 'hz=0');
            else
                lz = line(hx, [0 0], 'linestyle','-.','color','r','Tag', 'hz=0');
            end
        end
        
        if ~isempty(ht)
            hold on;
            set(ht,'Xdata', Xdata(TimeCurs), 'Ydata', Ydata(TimeCurs));
            hold off;
        end
        
        Ux = uicontextmenu;
        hux = xlabel(X.desp,'uicontextmenu',Ux);
        set(hux,'fontsize',FS(8));
        
        Uy = uicontextmenu;
        huy = ylabel(Y.desp,'uicontextmenu',Uy);
        set(huy,'fontsize',FS(8));
        
        xdesp = uimenu(Ux,'Label',X.help);
        ydesp = uimenu(Uy,'Label',Y.help);
        
        htt = title(titl,'uicontextmenu',Ux);
        set(htt,'fontsize',FS(8),'Tag','ttitle');
    else
        plot_res(h,X,Y,titl,D, 0);
    end
end

set(h,'userdata',D,'Tag','in_MD_plot','Visible','on','fontsize',FS(8));
