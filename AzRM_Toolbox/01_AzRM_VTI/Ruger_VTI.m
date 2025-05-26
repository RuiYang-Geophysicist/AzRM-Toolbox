function [obseis_Ruger] = Ruger_VTI(Tmax, thickness, dt, fdom, tlen, angle, rho, C_VTI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of synthetic seismic data in VTI media using Rüger's approximation
%
% INPUT:
%   Tmax        : Maximum modeling time (sec)
%   thickness   : Layer thicknesses (km)
%   dt          : Time sampling interval (sec)
%   fdom        : Dominant frequency of Ricker wavelet (Hz)
%   tlen        : Time length of the wavelet (sec)
%   angle       : Incidence angles (degrees)
%   rho         : Density (g/cm^3)
%   C_VTI       : Stiffness matrix of VTI model (6x6xn)
%
% OUTPUT:
%   obseis_Ruger: Synthetic seismic angle gather
%
% REFERENCE:
%   Rüger, A. (2002). Reflection coefficients and azimuthal AVO analysis.
%
%   Written by Rui Yang, Tongji University, 2024/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(thickness)
    C11(i,1)=C_VTI(1,1,i);C13(i,1)=C_VTI(1,3,i);C33(i,1)=C_VTI(3,3,i);
    C44(i,1)=C_VTI(4,4,i);C66(i,1)=C_VTI(6,6,i);
end

    vp=sqrt(C33./rho);
    vs=sqrt(C44./rho);
    
    epsilon=(C11-C33)./(2.*C33);
    gamma=(C66-C44)./(2.*C44);
    delta=((C13+C44).^2-(C33-C44).^2)./(2.*C33.*(C33-C44));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vp=1000*vp;                   % m/s
    vs=1000*vs;                   % m/s
    rho=1000*rho;                  % g/cm^3
    
    delta_T=0.002;
    nSamples = round(Tmax / delta_T) + 1;
    
    dz=1;
    zz=1000*[0;cumsum(thickness(1:end-1))];
    
    z=(0:dz:max(zz)+1000*thickness(end)*dz)';
    
    vp2=zeros(size(z));
    vs2=vp2;
    rho2=vp2;
    
    for k=1:length(zz)-1
        ind=find(z>=zz(k)&z<zz(k+1));
        vp2(ind)=vp(k);
        vs2(ind)=vs(k);
        rho2(ind)=rho(k);
        epsilon2(ind)=epsilon(k);
        delta2(ind)=delta(k);
    end
    
    ind=find(z>=zz(end));
    vp2(ind)=vp(end);
    vs2(ind)=vs(end);
    rho2(ind)=rho(end); 
    epsilon2(ind)=epsilon(end);
    delta2(ind)=delta(end);
    
    z2=z;
    dzvec=diff(z2);
    nz2=length(z2);
    loginttime=0;
    if(loginttime==0)
        %compute pp time curve
        tint=2*cumsum(dzvec./vp2(1:nz2-1));
    else
        %compute ps time curve
        tint=cumsum(dzvec./vp2(1:nz2-1)+dzvec./vs2(1:nz2-1));
    end

%Depth domain to time domain
    nt=round(tint(end)/delta_T)+1; % number of time samples needed
    vpt=zeros(nt,1);
    vst=zeros(nt,1);
    rhot=zeros(nt,1);
    epsilont=zeros(nt,1);
    deltat=zeros(nt,1); % 
    zt=zeros(nt,1);
    t1=0;
    
    for k=1:nt
        t2=t1+delta_T;
        ind=near(tint,t1,t2);
        vpt(k)=mean(vp2(ind));
        vst(k)=mean(vs2(ind));
        rhot(k)=mean(rho2(ind));
        epsilont(k)=mean(epsilon2(ind));
        deltat(k)=mean(delta2(ind));
        zt(k)=z2(ind(1));%use the first depth of each interval.
        t1=t2;
    end
    
    if(nt<=nSamples)
        npad=nSamples-nt;
        vpt(nt+1:nSamples)=vp(end);
        vst(nt+1:nSamples)=vs(end);
        rhot(nt+1:nSamples)=rho(end); 
        epsilont(nt+1:nSamples)=epsilon(end);
        deltat(nt+1:nSamples)=delta(end);
    end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ruger's approximation reflection coefficient of VTI medium
    k=zeros(nSamples,1);
    vp=vpt;
    vs=vst;
    den_pp=rhot;
    epsilon=epsilont;
    delta=deltat;
    
    k=vs.^2./vp.^2;
    
    nimpen=nSamples;
    nref=nimpen-1;
    
    [wavelet,~]=ricker(dt,fdom,tlen);
    
    [nw,nt]=size(wavelet);
    nw2=round(nw/2-1);
    mw=compute_w_matrix(wavelet,nref,nref,nw2);
 
    inc=angle./180*pi; % incident angle
    ninc=length(inc);
    
%Calculate the coefficient matrix G
    A0=zeros(nref,nref);
    B0=zeros(nref,nref);
    C0=zeros(nref,nref);
    D0=zeros(nref,nref);
    E0=zeros(nref,nref);
    Para_M0_w=[];
    for i=1:ninc
        for j=1:nref
            A0(j,j)=sec(inc(i))^2;
            B0(j,j)=-8*k(j)*sin(inc(i))^2;
            C0(j,j)=1-4*k(j)*sin(inc(i))^2;
            D0(j,j)=0.5*sin(inc(i))^2;
            E0(j,j)=0.5*sin(inc(i))^2*tan(inc(i))^2;
        end
        M0=[A0 B0 C0 D0 E0];
        M0_w=mw*M0;
        Para_M0_w=[Para_M0_w;M0_w];
    end
    G=Para_M0_w;
    
%calculate reflectivity
    Rvp=zeros(nref,1);
    Rvs=zeros(nref,1);
    Rden=zeros(nref,1);
    Rdelta=zeros(nref,1);
    Repsilon=zeros(nref,1);
    
    for i=1:nref
        Rvp(i)=(vp(i+1)-vp(i))/(vp(i+1)+vp(i));
        Rvs(i)=(vs(i+1)-vs(i))/(vs(i+1)+vs(i));
        Rden(i)=(den_pp(i+1)-den_pp(i))/(den_pp(i+1)+den_pp(i));
        Rdelta(i)=(delta(i+1)-delta(i));
        Repsilon(i)=(epsilon(i+1)-epsilon(i));
    end
    
    
    real_R=[Rvp;Rvs;Rden;Rdelta;Repsilon];
    trueseis=G*real_R;
    obseis_Ruger=reshape(trueseis,nref,ninc);
end


function [w,twave]=ricker(dtt,fm,tlength)
    if(nargin<3)
        tlength=127.*dtt;
    end
    if(nargin<2)
        fm=15.; 
    end
    n=round(tlength/dtt)+1;
    nzero=n/2;%zero time sample is here.
    nr=n-nzero-1;%number of samples to the right of nzero
    nl=nr;%number of samples to the right of nzero
    twave=dtt*(-nl:nr)';
    pf=pi^2*fm^2;
    w=(1-2.*pf*twave.^2).*exp(-pf*twave.^2);
end

function w1=compute_w_matrix(w,ns,nr,nw1)
    w1=zeros(ns,nr);
    nw=length(w);
    nw2=nw-nw1-1;
    
    for i=1:nr
        j1=nw1+2-i;
        w_max=nw2+i;
        j2=nw;
        if(j1<1)j1=1;end
        if(w_max>ns)
            w_max=ns;
            j2=w_max-i+1+nw1;
        end
        
        for j=j1:j2
            w1(i+j-nw1-1,i)=w(j);
        end
    end   
end

function ind=near(v,val1,val2)
 if nargin<3
    val2=val1;
 end
 
 %allow for the syntax near(val1,val2,v);
 if(length(val2)>1 && nargin>2)
    tmp=v;
    tmp2=val1;
    v=val2;
    val1=tmp;
    val2=tmp2;
  end
    
  ilive=find(~isnan(v));
  test=abs(v(ilive)-val1);
  L1= test==min(test);
  test=abs(v(ilive)-val2);
  L2= test==min(test);
  L1=ilive(L1);
  L2=ilive(L2);

 	if L1<=L2
  		ind=min(L1):max(L2);
	else
    	ind=max(L1):-1:min(L2);
 	end
end

