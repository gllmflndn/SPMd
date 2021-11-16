function test_spmd

GLM_PATH = getenv('GLM_PATH');

spm defaults fmri
spm_jobman initcfg

matlabbatch{1}.spm.tools.spmd.compute.spmmat = {fullfile(GLM_PATH,'SPM.mat')};
matlabbatch{1}.spm.tools.spmd.compute.outdir = {''};

spm_jobman('run',matlabbatch);
