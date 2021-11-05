%$Id: spmd_MDTS_plot.m,v 1.3 2005/07/25 19:43:39 huizhang Exp $
function spmd_MDTS_plot(h,type)
%------------------------------------------------------------------------
% Draw time series plots for model detail
% h{1} - figure number
% h{2} - userdata for plot
% type - which one the user specified
%        1- whitened residual
%        2- raw residual
%------------------------------------------------------------------------
global D;
ht = findobj(h,'Tag','TS_in_MD_plot');
ht = sort(ht);

if type==1
    % draw time series plots for whitened residual
    spmd_mtsview('update',ht,...
		   [D(4:6).data],{D(4:6).desp},...
	       [D(9).data D(3).data ],{D(9).desp D(3).desp});
    
else
    % draw time series plots for raw residual
    spmd_mtsview('update',ht,...
		   [D(11:13).data],{D(11:13).desp},...
	       [D(9).data D(10).data ],{D(9).desp D(10).desp});
    
end
return
