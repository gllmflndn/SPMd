%
% FORMAT spmd_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
% Individual users can make copies which can be stored in their own
% matlab subdirectories. If ~/matlab is ahead of the SnPM directory
% in the MATLABPATH, then the users own personal defaults are used.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id: spmd_defaults.m,v 1.3 2004/11/19 19:43:16 huizhang Exp $

global SPMd_defs

% Spatial extent threshold parameters
%------------------------------------------------------------------------
SPMd_defs.MaxClkDist = 0.03;  % Maximum distance away from a point in
                              % order for that point to be clicked
SPMd_defs.MarkerSize =    3;  % Using for marker size
SPMd_defs.Marker     =  'o';  % Using for marker style
