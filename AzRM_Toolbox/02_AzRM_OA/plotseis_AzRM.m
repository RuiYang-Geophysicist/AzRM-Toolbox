function fact=plotseis_AzRM(smat,t,x,vaflag,fact,clplvl,flipy,kolor,hax)
% PLOTSEIS: plot a seismic trace gather using WTVA format
%
% fact=plotseis(smat,t,x,vaflag,fact,clplvl,flipy,kolor,hax)
%
% PLOTSEIS does a quick plot of a seismic matrix in the current axis of the
% current figure window. It plots the columns as wiggle traces centered at
% their x coordinates.  Type plotseis with no arguments for a demo.
% Alternatives are PLOTSEISMIC which creates a new window and has UI
% controls for scaling, phase rotation, etc; nd PLOTIMAGE. If you have more
% than 100 traces, PLOTIMAGE is usually preferred.
%
% smat ... the seismic matrix to be plotted. Traces are assumed stored in
%          the columns smat.
% t   ...  time coordinates of traces.
%          default is 1:nrows where nrows = number of rows in smat
% x ...    x coordinates of the traces
%          default is 1:ncols where ncols=number of columns 
% vaflag   if 0, then traces are plotted as wt (wiggle traces) 
%            if 1, then the traces are plotted wtva (shaded peaks)
%            if 2 then reverve polarity wtva (shaded troughs)
%          default is 1
% fact ... scaling factors. Make fact(1) bigger for bigger wiggle
%          traces. fact(2) controls the overall plot scale. If not
%          provided, it is computed as max(abs(smat)). To scale two
%          plots the same, capture the return value from the first and
%          provide is as fact(2) for the second. See "scaling" below. The
%          second value of fact may also be set to one of two special negative
%          values that act as flags for special behavior. fact(2)=-1 is similar
%          defaulting fact(2) except that max(abs(smat(:))) is declared as a
%          global variable that can be used to control other plots. Setting
%          fact(2)=-2 causes a search to be made for a global variable for
%          scaling. If one is found, it is used and the resulting plot will be
%          in true scale relative to the plot that published the global. If none
%          is found, the max(abs(smat(:))) is used.
%          default is [1.5 max(abs(smat(:)))]
% clplvl   clip level specified as a number of trace spacings. Trace excursions
%	       larger than this many trace spacings get clipped.
%          default is 1 
% flipy    if 1, then the y axis is reversed so that it goes from top
%          of window to the bottom
%	       default is 1 
% kolor    color to plot the traces
%          default is [0 0 0] (black) 
% hax ...  handle of the axis to plot the seismic matrix in. If -1, then a 
%          new figure is created.
%          default = gca
%
% NOTE: To scale two plots with respect to the maximum absolute value
% on the first plot, capture the return value from the first
% and provide it as the fifth argument for the second plot:
%	fact=plotseis(smat1,t1,x1);
%	plotseis(smat2,t2,x2,1,fact);
%
% Scaling: First, the entire seismic matrix is divided by fact(2). Thus, if fact(2)
%    is defaulted, the matrix is normalized to have its maximum absolute value be
%    unity. Then, for n traces, the x axis range is divided into n bins and each trace
%    is plotted within its corresponding bin. Trace excursions outside cliplvl*binwidth 
%    boundaries are clipped. Let trc be the kth trace (i.e. the kth column of smat).
%    Prior to plotting trc in the kth bin it is scaled as trc=trc*binwidth*fact(1)
%    Thus each trace is plotted such that if an amplitude of
%    unity is encountered it causes an excursion of fact(1) binwidths.
%    So, the default for fact of [1.5 max(abs(smat(:)))] will cause the
%    maximum absolute value of smat to deflect 1.5 binwidths and so be 
%    clipped at 2/3 of its peak.
%
% G.F. Margrave, Aug 1995-2014
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE
global PLOTSEISFACT2

if(nargin<1)
	% do a demo
		%Make a fake reflectivity
		t=0:.002:1.0;
		r=randn(size(t)).^5;
		%make a ricker wavelet
		tw=-.1:.002:.1;
		arg=(pi*15*tw).^2;
		w=(1-2.*arg).*exp(-arg);
		%convolve
		s=conv(r,w);
		s=s(51:length(t)+50)';
		s=s/max(s); %normalize
		
		smat=s*ones(1,20);
end
if nargin < 9 || isempty(hax)
	hax=gca;
end
if nargin<8 || isempty(kolor)
	kolor=[0 0 0];
end
if nargin<7 || isempty(flipy)
	flipy=1;
end
if nargin<6 || isempty(clplvl)
	clplvl=1;
end
if nargin<5 || isempty(fact)
	fact=1.5;
end
if nargin<4 || isempty(vaflag)
	vaflag=1;
end
if nargin<3 || isempty(x)
	ncols=size(smat,2);
	x=1:ncols;
end
if nargin<2 || isempty(t)
	nrows=size(smat,1);
	t=1:nrows;
end

if(length(x)>1)
	bnds=(max(x)-min(x))/(length(x)+1);
else
	error('invalid x coordinate specification');
end

if( length(fact)<2 )
    % Assume a well behaved matrix free of +inf or -inf.  NaNs will be
    % ignored
    s=max(abs(smat(:)));
    if isinf(s)
        % Hmm.  Weed out the +inf and -inf values and try again.
        ilive=~isinf(smat);
        s=max(abs(smat(ilive)));
    end
	fact=[fact s];
end

if( fact(2) == -1)
    ilive=~isnan(smat);
	fact(2)=max(abs(smat(ilive)));
    PLOTSEISFACT2=fact(2);
end
if( fact(2) == -2)
    if(isempty(PLOTSEISFACT2))
        ilive=~isnan(smat);
	    fact(2)=max(abs(smat(ilive)));
    else
        fact(2)=PLOTSEISFACT2;
    end
end

if( hax==-1 )
	figure;
else
	hfig=get(hax,'parent');
	figure(hfig);
	set(hfig,'currentaxes',hax);
end

if(flipy)
	set(gca,'ydir','reverse');
end

for k=1:length(x)
	trc=smat(:,k)/fact(2); %normalize
	ilive=find(~isnan(trc));
	%m=mean(trc(ilive));
	%trc=(trc-m)*bnds*fact(1); %scale into x units
    trc=trc*bnds*fact(1); %scale into x units
    trc=clip(trc,clplvl*bnds)+x(k);%clip and shift
	if(~vaflag)
		line(trc(ilive),t(ilive),'color',kolor);
	elseif( vaflag==1)
		%replace NaN's with mean values
		if(length(ilive)< length(trc))
			idead=find(isnan(trc));
			trc(idead)=x(k)*ones(size(idead));
		end
		wtva(trc,t,kolor,x(k),1,1,1);
	elseif(vaflag==2)
		wtva(trc(ilive),t(ilive),kolor,x(k),-1,1,1);
	end
end

%simplezoom
%set(gcf,'name','Simple zooming installed, use MB1')
dx=abs(x(2)-x(1));
axis([min(x)-dx max(x)+dx min(t) max(t)])
%put the x axis on bottom
set(gca,'xaxislocation','bottom')

function [h,hva]=wtva(s,t,kolor,snot,polarity,dir,resampfact)
	if(nargin == 0)
	% do a demo
		%Make a fake reflectivity
		t=0:.002:2.0;
		r=randn(size(t)).^5;
		%make a ricker wavelet
		tw=-.1:.002:.1;
		arg=(pi*15*tw).^2;
		w=(1-2.*arg).*exp(-arg);
		%convolve
		s=conv(r,w);
		s=s(51:length(t)+50);
		s=s/max(s); %normalize
		%open up a figure
		figure
	end
	if(nargin<7)
		resampfact=1;
	end
	if(nargin<6)
		dir=1;
	end
	
	if(nargin<5)
		polarity=1;
	end

	if(nargin<4)
		ilive= ~isnan(s);
		snot=mean(s(ilive));
	end
	
	if(nargin<3)
		kolor='g';
	end
	
	if(abs(polarity)~=1)
		error('invalid polarity value');
	end
	
	%find the parts to shade
	% Method:
	%		1) Find all points not to be shaded
	%		2) find the boundary points of the groups of 1). That is
	%			determine how the points fall into contiguous groups and
	%			where these groups start and end
	%		3) Find any groups of 1) which are single isolated points and
	%			add an extra point so that the minimum group size is 2
	%		4) For all groups, interpolate in the zero crossings at the 
	%			beginning and end of the group and assign the first and
	%			last points in the group to these values
	%		5) assign all other points in the group to the value 0 (or snot)
	%		6) plot the polygon using patch with no edgecolor
	%

	%test for dead trace
	test=sum(abs(s-snot));
    %test for no peaks
%     itest=find(s>snot);
    
	if(test<1000*eps)
        ss=s;
		tt=t;
    else
		% resample
		if(resampfact>1)
			% my resampling
			[s,t]=resamp(s,t,(t(2)-t(1))/resampfact,[t(1) t(length(t))],0);
			% SIgnal toolbox resampling
			%s=resample(s,resampfact,1);
			%tmax=t(1)+(length(s)-1)*(t(2)-t(1))/resampfact;
			%t=t(1):(t(2)-t(1))/resampfact:tmax;
		end
		
	
		dt=t(2)-t(1);
		ss=polarity*[snot;s(:);snot];
		snot=polarity*snot;
		n=length(ss);
		tt=[t(1)-dt;t(:);t(length(t))+dt];
		% put snot on both ends. This forces the detection of isolated singular
		% points at the ends
		
		%find all points which won't be shaded
		il=find(ss<=snot);
		ilow=[0;il;n+1]; %end conditions
		%now find the beginning and ends of these zones
		ind=diff(ilow);
		ibdy=find(ind>1);
        if(~isempty(ibdy))
			% find singular groups
			ind= diff(ibdy)==1;
			ising=ilow(ibdy(ind)+1);
			if(~isempty(ising))
				if(ising(1)==0)
					ising(1)=[];
				end
				if(ising(length(ising))==n+1)
					ising(length(ising))=[];
				end
				if(~isempty(ising))
					for k=1:length(ising)
						%duplicate singular points
						ss=[ss(1:ising(k)); ss(ising(k)); ss(ising(k)+1:length(ss))];
						tt=[tt(1:ising(k)); tt(ising(k)); tt(ising(k)+1:length(tt))];
						ising=ising+1;
					end
				end
				% refind the points because the point count may have changed
				%find all points which won't be shaded
				il=find(ss<=snot);
				ilow=[0;il;n+1]; %end conditions
				%now find the beginning and ends of these zones
				ind=diff(ilow);
				ibdy=find(ind>1);
			end
        end
					
        if(~isempty(ibdy))
            ibdy1=ilow(ibdy);
            ibdy2=ilow(ibdy+1);
            if(ibdy1(1)==0)
                ibdy1(1)=[];
            end
            if(ibdy1(length(ibdy1))==n+1)
                ibdy1(length(ibdy1))=[];
            end
            if(ibdy2(1)==0)
                ibdy2(1)=[];
            end
            if(ibdy2(length(ibdy2))==n+1)
                ibdy2(length(ibdy2))=[];
            end
            
            %interpolate in zero crossings at the beginnings
            tnot=(snot-ss(ibdy1))./(ss(ibdy1+1)-ss(ibdy1));
            tnot=tnot.*(tt(ibdy1+1)-tt(ibdy1))+tt(ibdy1);
            
            tt(ibdy1)=tnot;
            %interpolate in zero crossings at the ends
            tnot=(snot-ss(ibdy2))./(ss(ibdy2-1)-ss(ibdy2));
            tnot=tnot.*(tt(ibdy2-1)-tt(ibdy2))+tt(ibdy2);
            tt(ibdy2)=tnot;
            
            % set all troughs to snot
            ss(il)=snot*ones(size(il));
            
            % make sure we have snot on both ends
            ss=[snot;ss(2:length(ss)-1);snot];
            tt=[t(1);tt(2:length(tt)-1);t(length(t))];
        else
            ss=snot*ones(size(s));
            tt=t;
        end
	end
	
	if(dir==1)
		hva=patch('xdata',polarity*ss,'ydata',tt,...
			'edgecolor','none','facecolor',kolor);
	else
		hva=patch('xdata',tt,'ydata',polarity*ss,...
		'edgecolor','none','facecolor',kolor);
	end

	%plot with line
	if(dir==1)
		h=line(s,t,'color',kolor);
	else
		h=line(t,s,'color',kolor);
	end

	if(nargout==1)
		h=[h hva];
		clear hva;
    end
end

function trout=clip(trin,amp);
% find the samples to be clipped
 indicies=find(abs(trin)>amp);
% clip them
 trout=trin;
 trout(indicies)=sign(trin(indicies))*amp; 
end

end








