% The first version is done by Wenlin
% $Id: spmd_SD.m,v 1.4 2007/03/16 14:51:24 nichols Exp $

function spmd_SD(varargin)
% spmd_SD: create series of residual images. If residual images are
% not available, series of raw images will be displayed.
%=======================================================================
% - FORMAT specifications for the embedded callback functions
%=======================================================================
% FORMAT spmd_SD
% Create spatial detail window, 'Tag'ged 'Temporal Series'. It displays
% mean images and five serial images of residuals or raw data around the
% temporal cursor. 
%________________________________________________________________________
% Warning:
% 1. To use this function, several files should exist in the current
% working directory: SPM.mat, and SPMres.mat.
% 2. Before opening the Temporal Series window, temporal summary  window
% should be opened.
%________________________________________________________________________
%
%
% %W% %E%

%------------------------ Functions Called -----------------------------
% spm_str_manip
% spmd_check_registration
% spmd_orthviews
%-----------------------------------------------------------------------

global TimeCurs
global Data
global STR
global RNG

TP = TimeCurs;

if isempty(TP),
  error('There is no temporal viewer exist!');
end

try
  load SPMd_SD
catch
  load SPM
  VY = SPM.xY.VY;
end

if isempty(Data);
  if exist('PstRes') == 1,
    Data = PstRes;
    STR = 'StdRes';
    RNG = [-3 3];
  else 
    nScan = length(VY);
    fnm = VY(1).fname;
    Path = spm_str_manip(fnm,'H');
    str1 = 'There is no studentized residual image available in the data directory:';
    str2 = ['  (' Path ')'];
    str3 = 'Do you want to specify the residual images?';
    str4 = '(If NO, images of raw data will be used in spatial detail!)'; 
    %spm('alert!',{str1,str2,str3,str4},'Warning');
    tmp = spm_input({str1,str2,str3 str4},1,'bd','Yes|No',[1,0]);
    close(gcf)
    if tmp == 1,
      PstRes = spm_select(nScan,'^e.*\.img','Select residual images...');
      Data = PstRes;
      STR = 'StdRes';
      RNG = [-3 3];
      save SPMres PstRes;
    elseif tmp == 0,
      Data = str2mat(VY.fname);
      STR = 'Raw Data ';
      RNG = {'auto','auto','auto','auto','auto'};
    end;
  end;
end;

nScan = length(Data);
	
if (TP<1), TP = 1; end
	
%=======================================================================
%-Create temporal series of images
%=======================================================================
if TP==1
  sptl('img',1,'Mean.img',...
  'img',2:6,cellstr(Data(TP:TP+4,:)),...
  'Wind',1,'auto',...
  'Wind',2:6,RNG,...
  'Desc',1,'Mean Image',...	  
  'Desc',2:6,strcat(STR,{num2str(TP),num2str(TP+1),num2str(TP+2),num2str(TP+3),num2str(TP+4)}))
	
elseif TP==2
  sptl('img',1,'Mean.img',...
  'img',2:6,cellstr(Data(TP-1:TP+3,:)),...
  'Wind',1,'auto',...
  'Wind',2:6,RNG,...
  'Desc',1,'Mean Image',...	  
  'Desc',2:6,strcat(STR,{num2str(TP-1),num2str(TP),num2str(TP+1),num2str(TP+2),num2str(TP+3)}))
	
elseif TP==nScan
  sptl('img',1,'Mean.img',...
  'img',2:6,cellstr(Data(TP-4:TP,:)),...
  'Wind',1,'auto',...
  'Wind',2:6,RNG,...
  'Desc',1,'Mean Image',...	  
  'Desc',2:6,strcat(STR,{num2str(TP-4),num2str(TP-3),num2str(TP-2),num2str(TP-1),num2str(TP)}))

elseif TP==nScan-1
  sptl('img',1,'Mean.img',...
  'img',2:6,cellstr(Data(TP-3:TP+1,:)),...
  'Wind',1,'auto',...
  'Wind',2:6,RNG,...
  'Desc',1,'Mean Image',...	  
  'Desc',2:6,strcat(STR,{num2str(TP-3),num2str(TP-2),num2str(TP-1),num2str(TP),num2str(TP+1)}))
	
else	 
  sptl('img',1,'Mean.img',...
  'img',2:6,cellstr(Data(TP-2:TP+2,:)),...
  'Wind',1,'auto',...
  'Wind',2:6,RNG,...
  'Desc',1,'Mean Image',...	  
  'Desc',2:6,strcat(STR,{num2str(TP-2),num2str(TP-1),num2str(TP),num2str(TP+1),num2str(TP+2)}))

end


function sptl(varargin)
%------------------------------------------------------------------------
% Create spatial detail window.
% Basically, this subfunction is the same as spmd_MS.m without default
% images to be shown.
%------------------------------------------------------------------------

%load SPM

Imgs = cell(1,24);
Desc = cell(1,24);
Wind = cell(1,24);

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

% fire up orthoviews  (w/ colorbars w/ windowing)
spmd_check_registration(str2mat(Imgs{:}),'Scan Detail');
fg = findobj('Tag','Scan Detail');
set(fg,'numbertitle','off');
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
  set(get(h,'YLabel'),'String',Desc{i});
end

return
