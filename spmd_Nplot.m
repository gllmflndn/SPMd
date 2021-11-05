function spmd_Nplot(y,h,D)
% create and update normal probability plot
% FORMAT spmd_Nplot(y,h,D)
%    Y- data
%    h- axes handle to plot the QQ-plot
%    D- all information attached to the plot
%________________________________________________________________________
% @(#)spmd_Nplot.m


global TimeCurs

if nargin < 2, h = gca; end

%-Check if input is a vector.
%------------------------------------------------------------------------
[m,n]=size(y);
if m~=1 && n~=1
    error('The input data must be a vector');
elseif m==1
    y=y';
    m=n;
end

%-Create the normal probability of the quantiles.
%------------------------------------------------------------------------
eprob = [0.5./m:1./m:(m - 0.5)./m];
x  = spm_invNcdf(eprob,0,1)';
[sy, si] = sort(y);
[stmp, osi] = sort(si);

%-Find the first and third quartiles of the data.
%------------------------------------------------------------------------
q1y = spm_percentile(sy,25);
q3y = spm_percentile(sy,75);
medy = median(sy);
q1x = spm_invNcdf(0.25,0,1);
q3x = spm_invNcdf(0.75,0,1);
medx = median(x);

%-Find the intercept and slope of the reference line in the QQ-plot.
%------------------------------------------------------------------------
slope = (q3y-q1y)./(q3x-q1x);
intcept = medy - slope.*medx;

%-Calculate the two extreme points of the reference line
%------------------------------------------------------------------------
my = [max(sy) min(sy)];
mx = [(max(sy)-intcept)./slope (min(sy)-intcept)./slope];

%-window Parameters
%-----------------------------------------------------------------------
WS     = spm('WinScale');
FS     = spm('FontSizes');
PF     = spm_platform('fonts');

%-Create the normality probability plot.
%-----------------------------------------------------------------------
axes(h);
%plot(x,sy,'.',mx,my,'r-.');
plot(x(osi),y,'.',mx,my,'r-.');
hold on;


ht = plot(x(osi(TimeCurs)), y(TimeCurs),'r*');
set(ht,'Tag','SelPt');
hold off;

Lx = xlabel('Expected quantile');
hx = uicontextmenu;
Ly = ylabel('Observed data');
hy = uicontextmenu;
Lx = xlabel('Expected quantile');
hx = uicontextmenu;
Ly = ylabel('Observed data');
hy = uicontextmenu;

set(Lx,'uicontextmenu',hx,'fontsize',FS(8));
set(Ly,'uicontextmenu',hy,'fontsize',FS(8));
xdes = uimenu(hx,'Label','Expected quantiles of normal distribution.');
ydes = uimenu(hy,'Label','Sorted observed data in ascending order.');

ht = title('Normal Probability Plot');
set(ht,'fontsize',FS(8));
set(h, 'fontsize',FS(8), 'userdata',D);
