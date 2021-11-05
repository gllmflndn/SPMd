function spmd_Nplot(y,h,D)
% Create and update normal probability plot
% FORMAT spmd_Nplot(y,h,D)
% Y    - data
% h    - axes handle to plot the QQ-plot
% D    - all information attached to the plot
%__________________________________________________________________________
% @(#)spmd_Nplot.m


global TimeCurs

if nargin < 2, h = gca; end

%-Check if input is a vector
%--------------------------------------------------------------------------
[m,n] = size(y);
if m~=1 && n~=1
    error('The input data must be a vector.');
elseif m == 1
    y = y';
    m = n;
end

%-Create the normal probability of the quantiles
%--------------------------------------------------------------------------
eprob    = 0.5/m : 1/m : (m-0.5)/m;
x        = spm_invNcdf(eprob,0,1)';
[sy, si] = sort(y);
[~, osi] = sort(si);

%-Find the first and third quartiles of the data
%--------------------------------------------------------------------------
q1y  = spm_percentile(sy,25);
q3y  = spm_percentile(sy,75);
medy = median(sy);
q1x  = spm_invNcdf(0.25,0,1);
q3x  = spm_invNcdf(0.75,0,1);
medx = median(x);

%-Find the intercept and slope of the reference line in the QQ-plot
%--------------------------------------------------------------------------
slope  = (q3y-q1y) ./ (q3x-q1x);
intcept = medy - slope.*medx;

%-Calculate the two extreme points of the reference line
%--------------------------------------------------------------------------
my = [max(sy) min(sy)];
mx = [(max(sy)-intcept)./slope (min(sy)-intcept)./slope];

%-Create the normality probability plot.
%--------------------------------------------------------------------------
plot(h,x(osi),y,'.',mx,my,'r-.');
hold(h,'on');
ht = plot(h,x(osi(TimeCurs)), y(TimeCurs),'r*');
set(ht,'Tag','SelPt');
hold(h,'off');

FS = spm('FontSizes');
Lx = xlabel(h,'Expected quantile');
hx = uicontextmenu;
Ly = ylabel(h,'Observed data');
hy = uicontextmenu;
set(Lx,'uicontextmenu',hx,'FontSize',FS(8));
set(Ly,'uicontextmenu',hy,'FontSize',FS(8));
uimenu(hx,'Label','Expected quantiles of normal distribution');
uimenu(hy,'Label','Sorted observed data in ascending order');

ht = title(h,'Normal Probability Plot', 'FontSize',FS(8));
if nargin > 2, set(h, 'UserData',D); end
