function cfg = tbx_cfg_spmd
% MATLABBATCH Configuration file for toolbox 'SPMd'

if ~isdeployed, addpath(fileparts(mfilename('fullpath'))); end

%--------------------------------------------------------------------------
% spmd_compute Compute SPMd
%--------------------------------------------------------------------------
spmd_compute      = cfg_exbranch;
spmd_compute.tag  = 'compute';
spmd_compute.name = 'Compute SPMd';
spmd_compute.val  = cfg_spmd_compute; % add @ for lazy-loading
spmd_compute.help = {
    'Compute SPMd'
    };
spmd_compute.prog = @run_spmd_compute;
spmd_compute.vout = @vout_spmd_compute;

%--------------------------------------------------------------------------
% spmd_visualise Visualise SPMd
%--------------------------------------------------------------------------
spmd_visualise      = cfg_exbranch;
spmd_visualise.tag  = 'visualise';
spmd_visualise.name = 'Visualise SPMd';
spmd_visualise.val  = cfg_spmd_visualise; % add @ for lazy-loading
spmd_visualise.help = {
    'Visualise SPMd'
    };
spmd_visualise.prog = @run_spmd_visualise;
spmd_visualise.vout = @vout_spmd_visualise;

%--------------------------------------------------------------------------
% spmd SPMd
%--------------------------------------------------------------------------
cfg        = cfg_choice;
cfg.tag    = 'spmd';
cfg.name   = 'SPMd';
cfg.help   = {
    'SPMd is a toolbox for SPM which you can use to establish the validity'
    'of inferences in fMRI modeling through diagnosis of linear model'
    'assumptions, and to characterize fMRI signal and artifacts through'
    'exploratory data analysis.'
    }';
cfg.values = {spmd_compute, spmd_visualise};


%==========================================================================
function varargout = cfg_spmd_compute

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

%--------------------------------------------------------------------------
% spmmat Select SPM.mat
%--------------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {
    'Select the SPM.mat file that contains the design specification.'
    }';
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

%--------------------------------------------------------------------------
% outdir Output Directory
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory';
outdir.val     = {{''}};
outdir.help    = {
    'Specify the output directory.'
    'If no directory is given, files will be written in the same directory than the input 4D file.'
    }';
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];

[cfg,varargout{1}] = deal({spmmat, outdir});


%==========================================================================
function varargout = cfg_spmd_visualise

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

%--------------------------------------------------------------------------
% spmdmat Select SPMd.mat
%--------------------------------------------------------------------------
spmdmat         = cfg_files;
spmdmat.tag     = 'spmdmat';
spmdmat.name    = 'Select SPMd.mat';
spmdmat.help    = {
    'Select the SPMd.mat file.'
    }';
spmdmat.filter  = 'mat';
spmdmat.ufilter = '^SPMd\.mat$';
spmdmat.num     = [1 1];

[cfg,varargout{1}] = deal({spmdmat});


%==========================================================================
function out = run_spmd_compute(job)

spm('FnBanner','spm_diagnosis.m');

%-Change to the analysis directory
%--------------------------------------------------------------------------
cwd = pwd;
try
    swd = spm_file(job.spmmat{1},'fpath');
    cd(swd);
    fprintf('%-40s: %30s\n','SPM directory',spm_file(swd,'short30'));   %-#
catch
    error('Failed to change directory %s.',swd)
end

%-Load SPM.mat file
%--------------------------------------------------------------------------
load(fullfile(swd,'SPM.mat'));
SPM.swd = swd;

%-Compute diagnostic statistics
%--------------------------------------------------------------------------
spmd_comp_MS(SPM);                   % Model summary - spatial
spmd_comp_SS;                        % Scan summary - temporal
spmd_MkResid(SPM,job.outdir{1});     % Scan detail - explore spatial detail

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#

%-Change back directory
%--------------------------------------------------------------------------
cd(cwd);

%-Output structure
%--------------------------------------------------------------------------
out.spmdmat = fullfile(swd,'SPMd.mat');


%==========================================================================
function out = run_spmd_visualise(job)

spm('FnBanner','spm_diagnosis.m');

%-Change to the analysis directory
%--------------------------------------------------------------------------
cwd = pwd;
try
    swd = spm_file(job.spmdmat{1},'fpath');
    cd(swd);
    fprintf('%-40s: %30s\n','SPMd directory',spm_file(swd,'short30'));  %-#
catch
    error('Failed to change directory %s.',swd)
end

%-Load SPM.mat file
%--------------------------------------------------------------------------
load(fullfile(swd,'SPMd.mat'));
SPMd.swd = swd;

%-Visualise diagnostic statistics
%--------------------------------------------------------------------------
%spmd_MS(SPMd);
%spmd_SS(SPMd);
%spmd_MD(SPMd);
%spmd_SD(SPMd);

%-Change back directory
%--------------------------------------------------------------------------
cd(cwd);

%-Output structure
%--------------------------------------------------------------------------
out = [];


%==========================================================================
function dep = vout_spmd_compute(job)
dep = [];


%==========================================================================
function dep = vout_spmd_visualise(job)
dep = [];
