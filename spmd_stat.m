% The first version is done by Wenlin
% $Id: spmd_stat.m,v 1.9 2008/02/20 12:19:02 nichols Exp $

function [stat,p,coef] = spmd_stat(action,res,X,para,df)
% calculate 
%  1. Dubin Watson statistic
%  2. Cook Weisberg score statistic
%  3. Shapiro-Wilk satatistic
%  4. Proportion of outlier
%  5. Cumulative periodogram statistic
% ______________________________________________________________________
% FORMAT spmd_stat(action,res,x,para)
%     action - option
%     res  - residual data
%     X    - design matrix
%     Para - parameters
%     df   - degree of freedom of error
% ____________________ Function called _________________________________
%
%     spm_detrend
%     spm_Xcdf
%     spm_invNcdf
%     spm_Ncdf
%     spm_Icdf
% ______________________________________________________________________       
% @(#)spmd_stat.m	1.4 Wen-Lin Luo 04/06/16

%load SPM
%df = SPM.xX.erdf;

[nScan nVox] = size(res);

if strmatch(action,'sw'),

else
    if size(X,1) ~= nScan,
      error(['The length of the residuals and design matrix must be the' ...
	     ' same!']);
    end
end

ResSS = sum(res.^2);
    
switch lower(action)
  
case 'dw'  
%----------------------------------------------------------------------------
% Durbin Watson statistics
%----------------------------------------------------------------------------
  if isempty(para)
    %--------------------------------------------------------------------------
    % Calculate the approximate distribution parameters of D
    %--------------------------------------------------------------------------
    nDWBeta = size(X,2);
    
    DWD1X = diff(X);
    DWD2X = diff(DWD1X);
    
    DWDP  = (2*nScan-1)-trace(DWD1X'*DWD1X*(inv(X'*X)));
    DWDQ  = 2*(3*nScan-4) - 2*trace(DWD2X'*DWD2X*(inv(X'*X))) + ...
	   trace((DWD1X'*DWD1X*(inv(X'*X)))^2);    
    DWED  = DWDP/(nScan-(nDWBeta)-1);
    DWVarD= 2/((nScan-(nDWBeta)-1)*(nScan-(nDWBeta)+1))*(DWDQ-DWDP*DWED);
    
    DWp   = (DWED*(4-DWED)/DWVarD-1)*DWED/4;
    DWq   = (DWED*(4-DWED)/DWVarD)-1-DWp;
    
    coef  = [DWp;DWq];
  else
    coef  = para;
  end
  
    %-----------------------------------------------------------------
    %-Durbin-Watson statistic
    %-----------------------------------------------------------------
    DW   = sum((diff(res)).^2)./ResSS;	%-Durbin-Watson statistic
    stat = coef(1)*(4-DW)./(coef(2)*DW);
    p	 = 1-betainc(1-DW/4,coef(2), coef(1));	%-DW probability
    p	 = p+(p==0)*eps;
    p	 = -log10(p);
    stat = {stat};
case 'score'    
%-----------------------------------------------------------------------
% Score test / Cook Weisberg score statistic
%-------------------------------------------------------------------------
  SSr   = nScan*(res.^2)*(spdiags((ResSS.^(-1))',0,nVox,nVox)); 
  
  meanSSr = mean(SSr);
  [nrowSSr ncolSSr] = size(meanSSr);
  
  if ncolSSr > 1
      SYY	= sum((SSr-repmat(meanSSr,nScan,1)).^2); %-scaled squared residual
  else
      SYY	= sum((SSr-repmat(meanSSr(1,1),nScan,1)).^2); %-scaled squared residual
  end
  

  coef  = para;
    %-----------------------------------------------------------------
    % regress the scaled squared residual on possible regressors
    %-----------------------------------------------------------------
    if isempty(coef),
      Rres  = SSr-X*pinv(X)*SSr;
      %def   = size(X,2)-1;
      def   = max(1,size(X,2)-1); %add Tom's change
    elseif coef == 1,
      X     = spm_detrend(X);
      SSr   = spm_detrend(SSr);
      
      sumX2 = sum(X.^2)'; % define
      [xindex, yindex] = size(sumX2);
      for i=1:xindex
          for j=1:yindex
              if abs(sumX2(i,j)-0)<1e-7, sumX2(i,j)=0.000001; end
          end
      end
              
      
      %Rres  = SSr - X*spdiags(sum(X.*SSr)'./sum(X.^2)', ...
      %0,nVox,nVox);
      Rres = SSr - X*spdiags(sum(X.*SSr)'./sumX2,0,nVox,nVox);
                  
      def   = 1;
    end
    
    stat = (SYY-sum(Rres.^2))./2;		%-Score statistics
    p	 = 1-spm_Xcdf(stat,def); 	        %-Score probability  
    p	 = p+(p==0)*eps;
    p	 = -log10(p);
    stat = {stat};
    
case 'sw'  
%-------------------------------------------------------------------------
%-For Shapiro-Wilks Normality Test
%-------------------------------------------------------------------------
  if isempty(para)
    %-----------------------------------------------------------------
    % Calculate necessary coefficient for Shapiro-Wilks test
    %-----------------------------------------------------------------  
      
      % Calculate expected normal order statistics
      %-----------------------------------------------------------------  
      eprob = [0.5./nScan:1./nScan:(nScan-0.5)./nScan];
      enorm = spm_invNcdf(eprob,0,1);
      
      % Calculate coefficient for W statistic
      %-----------------------------------------------------------------     
      one = [0 ones(1,nScan-2) 0];
      SWa = 2*enorm;
      SWa1s = (6.*nScan+7)./(6.*nScan+13).*...
	      ((exp(1)./(nScan+2).*((nScan+1)./(nScan+2)).^(nScan-2)).^0.5);
      SWa(:,1) = -((SWa1s./(1-2.*SWa1s))*(one*(SWa.^2)')).^0.5;
      SWa(:,nScan) = ((SWa1s./(1-2.*SWa1s))*(one*(SWa.^2)')).^0.5;
      SWa = SWa./(SWa*SWa').^0.5;
      
      %
      % Calculate some constants for Approximation for Significance level
      % of W: lambda, mean and SD for the transformation of W
      % poly = log(nScan)-5;     % for sample size greater than 20, d=5 
      %-----------------------------------------------------------------     
      SWcoeff1 = [0.002989646 0.00879701 -0.0241665 0 0.318828 0.480385];
      SWcoeff2 = [-0.01504614 -0.03513666 0.1066339 -0.04183209...
		  -1.37888 -1.91487];
      SWcoeff3 = [0.003852646 -0.03215018 -0.01638782 0.1773538...
		  -0.331885 -1.015807 -3.73538];
      SWlambda = polyval(SWcoeff1,log(nScan)-5);
      SWmu = exp(polyval(SWcoeff2,log(nScan)-5));
      SWsd = exp(polyval(SWcoeff3,log(nScan)-5));
      
      coef = {SWa SWlambda SWmu SWsd};
  else
    coef = para;
  end
  
    %-----------------------------------------------------------------
    %-Shapiro-Wilks statistics
    %-----------------------------------------------------------------
    sRes  = sort(res);
    WRes  = (coef{1}*sRes).^2./sum((res-repmat(mean(res),nScan,1)).^2);
    
    TWRes = (1-WRes).^coef{2};	   %-approximate normalizing
                                   % transformation of W 
    stat  = (TWRes-coef{3})./coef{4};	   %-standard normal deviate for
					   % transformation of W 
    p     = 1-spm_Ncdf(stat,0,1);  %-Shapiro-Wilks probability
    p     = p+(p==0)*eps;
    p     = -log10(p);
    stat = {stat};
 
case 'outl'    
%-----------------------------------------------------------------
%-outlier count
%-----------------------------------------------------------------
  if isempty(para),
    coef = 3;
  else
    coef = para;
  end
  
  MSE	= ResSS./df;		                %-Mean squared error    
  Outh	= X*pinv(X);		                %-leverage 
  Outh	= diag(Outh);		                %
  I_null= (1-Outh)<sqrt(eps);                   %-Zero variance indicates
                                                % null (exact zero) outliers
  Outh(I_null)=0;                               %-Residual zero, so this
                                                % can be whatever
  Sres	= res.*((1-Outh)*MSE).^(-0.5);    %-standardized residual
  Iout  = (abs(Sres) > coef);	                %-index of outliers
  
  %-spatial outlier count
  %---------------------------------------------------------------
  p     = 1-spm_Icdf(sum(Iout),nScan,2*(1-spm_Ncdf(coef)));
  p     = p+(p==0)*eps;
  p     = -log10(p);
  
  %-temporal outlier count
  stat = {sum(Iout) sum(Iout,2)};
case 'cp'    
%-----------------------------------------------------------------
%-cumulative periodogram
%-----------------------------------------------------------------
  nVar = rank(X);

  if isempty(para),
    %--------------------------------------------------------------
    %-Calculate the vector for uncorrelating the residuals
    %--------------------------------------------------------------
    if nVar == 0,
      coef  = 0;
    else
      Xb0   = X(1:nVar,:);
      Xb1   = X(nVar+1:nScan,:);
      
      [V,D] = eig(Xb0*pinv(X'*X)*Xb0');
      D     = spdiags(D,0);
      V     = V./(repmat(sqrt(sum(V.*V)),nVar,1));
      H     = sum(D==1);
      
      blus1 = zeros(nVar,nVar);
      for i = H+1:nVar
	blus1 = blus1 + D(i)/(1+D(i))*V(:,i)*V(:,i)';
      end
      
      coef = Xb1*pinv(Xb0)*blus1;
    end
  else
    coef = para;
  end

    %-Cumulative periodogram
    %-----------------------------------------------------------------
    power		 = zeros(floor((nScan-nVar-1)/2),nVox);
    
    for i = 1:nVox
      %- Calculate BLUS residual vector
      %---------------------------------------------------------------
      Rblus      = res(nVar+1:nScan,i)-coef*res(1:nVar,i);
      
      %- Calculate periodogram
      %---------------------------------------------------------------
      Spec       = abs(fft(Rblus,nScan-nVar)).^2;
      
      % Select first half  
      if rem(nScan-nVar,2),         % nfft odd
	select   = (2:(nScan-nVar+1)/2)';  %-don't include DC or Nyquist
                                           % components 
      else
	select   = (2:(nScan-nVar)/2)';
      end
      freq       = (select-1)*2/(nScan-nVar);
      
      % Calculate the single-sided spectrum which includes the full power
      power(:,i) = [2*Spec(select)];
    end

    %-Calculate cumulative periodogram
    CPseries    = cumsum(power);
    Flength     = size(CPseries,1);
    CPS	= CPseries./(repmat(sum(power),Flength,1));
    
    %-setup the upper and lower statistic for the Cumulative periodogram
    MaxA1 = abs(CPS(1:Flength-1,:)-...
		(repmat([1:Flength-1]',1,nVox))/(Flength));
    MaxB1 = abs(CPS(1:Flength-1,:)-...
		(repmat([1:Flength-1]',1,nVox)-1)/(Flength));
    
    %-Calculate the Cumulative periodogram test statistic
    stat  = max(max(MaxA1,MaxB1));
    
    %-Calculate the coresponding p-value for the CPD
    CPZ  = stat*(sqrt(Flength-1)+0.12+0.11/sqrt(Flength-1));
    p    = 1-sqrt(2*pi)*(sum(exp(-((2*repmat([1:30]',1,nVox)...
		  -1)*pi).^2./(8*repmat(CPZ,size([1:30],2),1).^2))))./CPZ;
    p    = p+(p==0)*eps;
    p    = -log10(p);
           
    %-Calculate sum of periodogram
    stat = {stat sum(power,2) freq};
 
otherwise,
 warning('Unknown action string');
end;

