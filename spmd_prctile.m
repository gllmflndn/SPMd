% The first version is done by Wenlin
% $Id: spmd_prctile.m,v 1.2 2004/12/27 16:34:37 huizhang Exp $


function q = spmd_prctile(x,p)
% spmd_prctile(x,p): calculate the p percentile of x.
% FORMAT q=spmd_prctile(x,p)
% Input:
% x 	- Either a matrix or a vector.
% p 	- Vector of percentiles (in [0 100])
% Output:
% q	- values of percentiles.
%_________________________________________________________________________
% @(#)spmd_prctile.m	1.3 Wen-Lin Luo 02/10/03

[pm pn]=size(p);
if ((pm ~= 1) & (pn ~= 1)),
  error('p must be a scalar or a vector');
end

if (min(p)<0 | max(p)>100),
  error('The percentile should be in [0 100]');
end;

[m n] = size(x);

if m == 1,
  if n == 1,
    q = x*ones(length(p),1);
    return
  else
    x=x';
  end
end

x  = sort(x);
L  = length(x);
pp = [0 100*(0.5:L-0.5)/L 100];
x  = [min(x);x;max(x)];

q = interp1(pp,x,p);

