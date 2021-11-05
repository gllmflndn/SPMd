function Ind = spmd_pointer(varargin)
% Function used to point out the clicked point.
% FORMAT I = spmd_pointer(F)
% I    Index of selected point
% F    Vector of axes handles to look for clicked point
%
% FORMAT spmd_pointer('Init',)
% Initializes pointer structures
%
% FORMAT spmd_pointer('Set',I)
% update pointer structures
%__________________________________________________________________________
% @(#)spmd_pointer.m


global TimeCurs
global SPMd_defs

set(gcf,'DoubleBuffer', 'on'); %- elimate flash

if isempty(SPMd_defs), spmd_defaults; end
MaxClkDist = SPMd_defs.MaxClkDist;

if nargin==0; end                                    %- do nothing

TR = varargin{1};
if (isempty(TR)) return;  end                        %- do nothing
Ind = [];

if ~isstr(TR)
    F = TR;
else
    switch(lower(TR))
        case 'init'
            %- create plot
            %--------------------------------------------------------------
            FC       = get(gcf,'children');
            A        = strcmp(get(FC,'type'),'axes');
            F        = findobj(FC(find(A)),'Tag','in_MD_plot');
            NumChild = length(F);
            
            for Index=1:NumChild
                axes(F(Index));
                hp = findobj(F(Index), 'Tag','SelPt');
                Cdata = get(F(Index),'children');
                Xdata = get(Cdata(end),'xdata');
                Ydata = get(Cdata(end),'ydata');
                OrigXdata = get(hp,'Xdata');
                OrigYdata = get(hp,'Ydata');
                
                if (~isempty(hp))
                    %-  SelPt is already there
                    if (OrigXdata ~= Xdata(TimeCurs)) || (OrigYdata ~= Ydata(TimeCurs))
                        %- no change in index
                        set(hp,'visible','on');
                        set(hp, 'Xdata', Xdata(TimeCurs), 'Ydata', Ydata(TimeCurs), 'Tag', 'SelPt');
                    end
                else
                    %- No select point
                    %- Add a new point
                    hold on
                    hp = plot(Xdata(TimeCurs),Ydata(TimeCurs),'r*');
                    set(hp,'Tag','SelPt');
                    hold off
                end
            end
            
            return
        case 'set'
            %- update plot
            %-------------------------------------------------------------
            Ind      = varargin{2};
            
            %- find 'in_MD_plot' plot
            FC       = get(findobj('Tag','SPMd_MD'),'children');
            A        = strcmp(get(FC,'type'),'axes');
            F        = findobj(FC(find(A)),'Tag','in_MD_plot');
            NumChild = length(F); %- number of such plots
            
            for Index=1:NumChild
                axes(F(Index));
                hp = findobj(F(Index), 'Tag','SelPt');
                Cdata = get(F(Index),'children');
                Xdata = get(Cdata(end),'xdata');
                Ydata = get(Cdata(end),'ydata');
                
                OrigXdata = get(hp,'Xdata');
                OrigYdata = get(hp,'Ydata');
                if (OrigXdata ~= Xdata(Ind)) || (OrigYdata ~= Ydata(Ind))
                    hold on;
                    set(hp, 'Xdata', Xdata(Ind), 'Ydata', Ydata(Ind), 'Tag', 'SelPt');
                    set(hp, 'visible','on');
                    hold off;
                end
            end
            return
    end
end


%- if the currentpoint is not within the range of data, then return
%-------------------------------------------------------------------
fg = gca;
p = get(fg,'currentpoint');

if (p(1,3) ~= -p(2,3))
    Ind = [];
    return
end


%- Find the (x,y) position for the current data point
%------------------------------------------------------------------
P = p(1,1:2);
Cdata = get(fg,'children');
Ctype = get(Cdata(end),'type');

if strcmp(Ctype,'text')
    Id = [];
    return;
else
    
    load SPM;
    nScan = SPM.nscan;
    
    Xdata = get(Cdata(end),'xdata');
    Ydata = get(Cdata(end),'ydata');
    
    [nrow, ncol] = size(Xdata);
    
    Index = 1;
    while (ncol~= nScan)
        Xdata = get(Cdata(end-Index),'xdata');
        Ydata = get(Cdata(end-Index),'ydata');
        [nrow, ncol] = size(Xdata);
        Index = Index + 1;
    end
    
    rangeX = get(fg,'xlim');
    scaleX = (rangeX(2) - rangeX(1))*10;
    
    rangeY = get(fg,'ylim');
    scaleY = (rangeY(2) - rangeY(1))*10;
    
    c = sqrt( (((Xdata-p(1,1))/scaleX).^2 + ((Ydata-p(1,2))/scaleY).^2)*100 );
    if min(c) >MaxClkDist
        Ind = [];
    else
        Ind = min(find(c == min(c)));
    end
    
    %- Put the text of position of current data point in all the subplots in
    %  current figure
    %-----------------------------------------------------------------------
    for i = 1:length(F)
        C = get(F(i),'children');
        hp = findobj(F(i), 'Tag','SelPt');
        
        if (~isempty(Ind))
            xdata = get(C(end),'xdata');
            ydata = get(C(end),'ydata');
            
            %-window Parameters
            %-------------------------------------------------------------------
            WS     = spm('WinScale');
            FS     = spm('FontSizes');
            PF     = spm_platform('fonts');
            set(hp,'Xdata', xdata(Ind), 'Ydata', ydata(Ind), 'Tag', 'SelPt');
            set(hp,'visible', 'on');
        else
            set(hp,'visible', 'off');
        end
    end
end
