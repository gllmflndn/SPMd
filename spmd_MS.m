% The first version is done by Wenlin
% $Id: spmd_MS.m,v 1.4 2005/03/25 17:22:02 huizhang Exp $

function spmd_MS(varargin)
%
% Creates page of spatial summary images
%=======================================================================
% - FORMAT specifications for the embedded callback functions
%=======================================================================
% FORMAT spmd_MS
% Six default images are shown in the spatial summary window.
%   F-statistic
%   -log10(P-value) of Durbin-Watson statistic
%   -log10(P-value) of Cook-Weisberg score statistic
%   -log10(P-value) of Shapiro-Wilk statistic
%   Mean image
%   Standard Deviation image
%
% FORMAT spmd_MS('Img',n,'FileNm')
% Specify the filenames of the images.
% n            - indices of the images.
% FileNm       - filenames of the images.
%
% FORMAT spmd_MS('Desc',n,'Image name')
% Specify brief description of the images.
% n            - indices of the images.
% Image name   - descriptions of the images.
%
% FORMAT spmd_MS('window',n,win)
% Window the intensity of the images
% n            - indices of the images.
% win          - a 1x2 vector for the range of intensities of the
%                images. 
%_________________________________________________________________________
% spmd_MS opens the spatial summary window of multiple statistical
% images. There are at most 24 images can be shown in the window. 
%     
%                         ---------------------
%
% Warning:
% The SPM.mat file should exist in the current working directory.
%     
%                         ---------------------
%
% Usage of this function:
% Multiple commands can be issued, like
%   spmd_MS('img',2,'spmF_0001','Window',2,[-4 4])  
%
% Multiple files can be specified like
%   spmd_MS('Img',1:3,{'Img1','Img2','Img3'})
%
% Multiple descriptions can be specified like
%   spmd_MS('Desc',1:3,{'Title1','Title2','Title3'})
%
% Multiple windows can be specified like
%   spmd_MS('Wind',1:3,[0 50])
% or
%   spmd_MS('Wind',1:3,{[0 50],[-6 6],[0 10]})
%________________________________________________________________________
% @(#)spmd_sptl.m	1.15 Tom Nichols & Wen-Lin Luo 03/07/16
%
% @(#)spmd_MS.m	1.2 04/07/08 

%------------------------ Functions Called -----------------------------
% spm_str_manip
% spmd_check_registration
% spmd_orthviews
%-----------------------------------------------------------------------

%load SPM

%========================================================================
%- Set up the default images
%========================================================================

nImgDef = 6;
% Default images

DefImgs = {'spmF_0002','PCorr',...
	'PHomo1','PNorm',...
	'Mean','ResRMS_MDC'};
% Default descriptions
DefDesc = {'F-statistics','Durbin-Watson (-log10 p)',...
	'Score (-log10 p)','Shapiro-Wilks (-log10 p)',...
	'Mean','Standard Deviation'};
% Default windows
DefWind = {'auto',[0 10],...
	[0,10],[0 10],...
	'auto',[0 1]};

%-set the maximal number of images shown in the spatial summary window.
%-----------------------------------------------------------------------
Imgs = cell(1,24);
Desc = cell(1,24);
Wind = cell(1,24);

%========================================================================
%- Set up different formats
%========================================================================
a = 1;
while (a <= nargin)
  
  switch lower(varargin{a})

   case 'img'
    if (nargin<3), error('Bad useage (did you include image number?)');end
    ks   = varargin{a+1};
    fNm = varargin{a+2};
    if (length(ks)>1 & (~iscell(fNm) | length(ks)~=length(fNm)))
      error('Must specify equal number of img nums and img names'); end
    for  k=1:length(ks)
      if (length(ks)==1)
	Imgs{ks} = fNm;
      else
	Imgs{ks(k)} = fNm{k};
      end
    end
    a = a + 3;

   case 'desc'
    if (nargin<3), error('Bad useage (did you include image number?)');end
    ks   = varargin{a+1};
    Dsc = varargin{a+2};
    if (length(ks)>1 & (~iscell(Dsc) | length(ks)~=length(Dsc)))
      error('Must specify equal number of img nums and descriptions'); end
    for  k=1:length(ks)
      if (length(ks)==1)
	Desc{ks} = Dsc;
      else
	Desc{ks(k)} = Dsc{k};
      end
    end
    a = a + 3;

   case {'window', 'wind', 'win'}
    if (nargin<3), error('Bad useage (did you include image number?)');end
    ks   = varargin{a+1};
    Win = varargin{a+2};
    if (~iscell(Win) & length(Win)==2)
      % Allow multiple images to be set to a common window
      Wind(ks) = {Win};
    else
      if (length(ks)>1 & (~iscell(Win) | length(ks)~=length(Win)))
	error('Must specify equal number of img nums and windows'); end
      for  k=1:length(ks)
	if (length(ks)==1)
	  Wind{ks} = Win;
	else
	  Wind{ks(k)} = Win{k};
	end
      end
    end
    a = a + 3;

   otherwise
    error('Unknown command');
    
  end

end

%========================================================================
%-set up the images to be displayed.
%========================================================================
if (nargin==0)
  % Set defaults
  for i=1:nImgDef
    Imgs{i} = DefImgs{i};
    Desc{i} = DefDesc{i};
    Wind{i} = DefWind{i};
  end
end

% Weed out unused elements
Gd = [];
for i=1:24
  if ~isempty(Imgs{i}), 
    Gd = [Gd i];
    if isempty(Desc{i}),
      Desc{i} = spm_str_manip(Imgs{i},'ts');
    end
    if isempty(Wind{i}),
      Wind{i} = 'auto';
    end
  end
end
Imgs = Imgs(Gd);
Desc = Desc(Gd);
Wind = Wind(Gd);

nImg = length(Imgs);

%========================================================================
%-Create the model summary window
%========================================================================
% fire up orthoviews  (w/ colorbars w/ windowing)
spmd_check_registration(str2mat(Imgs{:}),'Model Summary');
set(gcf,'DoubleBuffer','on');
fg = findobj('Tag','Model Summary');
set(fg,'numbertitle','off',...
    'name','Model Summary');
spmd_orthviews('fig',fg,'addcolorbar',1:nImg);
drawnow
for i=1:nImg
  if ~strcmp(Wind{i},'auto')
    spmd_orthviews('fig',fg,'window',i,Wind{i});
  end
end

st = spmd_orthviews('fig',fg,'getst');

% Name the images
for i=1:nImg
  h = st.vols{i}.ax{1}.ax;
  set(h, 'YAxisLocation','right');
  set(get(h,'YLabel'),'Interpreter','tex','String',Desc{i});
end


