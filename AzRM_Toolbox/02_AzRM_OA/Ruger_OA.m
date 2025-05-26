function [obseis_OA] = Ruger_OA(Tmax, thickness, dt, fdom, tlen, angle, model, phi1, phi2, delta_phi, phi0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute synthetic seismic gather in orthorhombic media using Ruger and Chen approximations
%
% INPUT:
%   Tmax        : Maximum modeling time (s)
%   thickness   : Layer thicknesses (km)
%   dt          : Time sampling interval (s)
%   fdom        : Dominant frequency of Ricker wavelet (Hz)
%   tlen        : Time length of wavelet (s)
%   angle       : Array of incidence angles (deg)
%   model       : Elastic model with columns: [thick, vp, vs, rho, epsilon, delta, gamma, e]
%   phi1, phi2  : Azimuth range (deg)
%   delta_phi   : Azimuth increment (deg)
%   phi0        : Fracture orientation (deg)
%
% OUTPUT:
%   obseis_OA   : Synthetic seismic gather [samples × azimuth]
%
% REFERENCES:
%   Rüger, A. (2002). Reflection coefficients and azimuthal AVO analysis.
%   Chen et al., GJI (2017). AVAZ inversion using elastic impedance and MCMC.
%
%   Written by Rui Yang, Tongji University, 2024/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------- Extract elastic parameters ----------------------
vp       = 1000 * model(:, 2);     % m/s
vs       = 1000 * model(:, 3);     % m/s
rho      = 1000 * model(:, 4);     % kg/m3
epsilon  = model(:, 5);
delta    = model(:, 6);
gamma    = model(:, 7);
e        = model(:, 8);           % Fracture weakness

g        = (vs ./ vp).^2;
dN       = 4 * e ./ (3 * g .* (1 - g));
dT       = 16 * e ./ (3 * (3 - 2 * g));
Ip       = rho .* vp;
Is       = rho .* vs;

% -------------------- Depth-to-Time Transform --------------------------
delta_T  = 0.002;
nSamples = round(Tmax / delta_T) + 1;

% Rotate azimuthal domain to fracture-aligned coordinates
phi1 = phi1 - phi0;
phi2 = phi2 - phi0;
nazimuth = round((phi2 - phi1) / delta_phi) + 1;
nangle   = length(angle);

dz  = 1;
zz  = 1000 * [0; cumsum(thickness(1:end-1))];
z   = (0:dz:(max(zz) + 1000 * thickness(end)))';

vp2 = zeros(size(z));
Ip2 = vp2; Is2 = vp2; rho2 = vp2;
epsilon2 = vp2; delta2 = vp2; dN2 = vp2; dT2 = vp2;

for k = 1:length(zz)-1
    ind = z >= zz(k) & z < zz(k+1);
    vp2(ind)       = vp(k);
    Ip2(ind)       = Ip(k);
    Is2(ind)       = Is(k);
    rho2(ind)      = rho(k);
    epsilon2(ind)  = epsilon(k);
    delta2(ind)    = delta(k);
    dN2(ind)       = dN(k);
    dT2(ind)       = dT(k);
end

% Extend last layer
ind = z >= zz(end);
vp2(ind)      = vp(end);
Ip2(ind)      = Ip(end);
Is2(ind)      = Is(end);
rho2(ind)     = rho(end);
epsilon2(ind) = epsilon(end);
delta2(ind)   = delta(end);
dN2(ind)      = dN(end);
dT2(ind)      = dT(end);

% Compute vertical two-way travel time curve
dzvec = diff(z);
tint  = 2 * cumsum(dzvec ./ vp2(1:end-1));
nt    = round(tint(end) / delta_T) + 1;

% Time sampling of elastic logs
Ipt = zeros(nt,1); Ist = Ipt; rhot = Ipt;
epsilont = Ipt; deltat = Ipt; dNt = Ipt; dTt = Ipt;
zt = Ipt;

t1 = 0;
for k = 1:nt
    t2 = t1 + delta_T;
    ind = near(tint, t1, t2);
    Ipt(k)      = mean(Ip2(ind));
    Ist(k)      = mean(Is2(ind));
    rhot(k)     = mean(rho2(ind));
    epsilont(k) = mean(epsilon2(ind));
    deltat(k)   = mean(delta2(ind));
    dNt(k)      = mean(dN2(ind));
    dTt(k)      = mean(dT2(ind));
    zt(k)       = z(ind(1));
    t1 = t2;
end

% Padding
if nt < nSamples
    Ipt(nt+1:nSamples)      = Ip(end);
    Ist(nt+1:nSamples)      = Is(end);
    rhot(nt+1:nSamples)     = rho(end);
    epsilont(nt+1:nSamples) = epsilon(end);
    deltat(nt+1:nSamples)   = delta(end);
    dNt(nt+1:nSamples)      = dN(end);
    dTt(nt+1:nSamples)      = dT(end);
end

% -------------------- Build Chen's Reflectivity Model --------------------------
k   = Ist.^2 ./ Ipt.^2;
nref = nSamples - 1;

[wavelet, ~] = ricker(dt, fdom, tlen);
nw2 = round(length(wavelet)/2 - 1);
mw  = compute_w_matrix(wavelet, nref, nref, nw2);

% Incidence and azimuth reshape
inc_org = angle * pi / 180;
ninc = length(inc_org);
azimuth = phi1:delta_phi:phi2;
phi_org = azimuth * pi / 180;
nazi = length(phi_org);
ntr  = ninc * nazi;

inc = repmat(inc_org, 1, nazi);
phi = repelem(phi_org, ninc);

% Construct matrix G (Eq. 8 in Chen et al., 2017)
A0 = zeros(nref,nref); B0 = A0; C0 = A0;
D0 = A0; E0 = A0; F0 = A0; G0 = A0;
Para_M0_w = [];

for i = 1:ntr
    for j = 1:nref
        A0(j,j) = sec(inc(i))^2;
        B0(j,j) = -8 * k(j) * sin(inc(i))^2;
        C0(j,j) = 4 * k(j) * sin(inc(i))^2 - tan(inc(i))^2;
        D0(j,j) = sin(inc(i))^2 * tan(inc(i))^2;
        E0(j,j) = sin(inc(i))^2;
        F0(j,j) = -(1 - 2 * k(j) * (sin(inc(i))^2 * sin(phi(i))^2 + cos(inc(i))^2))^2 ...
                  / (2 * cos(inc(i))^2);
        G0(j,j) = 2 * k(j) * sin(inc(i))^2 * cos(phi(i))^2;
    end
    M0    = [A0 B0 C0 D0 E0 F0 G0];
    M0_w  = mw * M0;
    Para_M0_w = [Para_M0_w; M0_w];
end

G = Para_M0_w;

% Compute reflectivity vector
real_R = [
    diff(Ipt) ./ (Ipt(1:end-1) + Ipt(2:end));
    diff(Ist) ./ (Ist(1:end-1) + Ist(2:end));
    diff(rhot) ./ (rhot(1:end-1) + rhot(2:end));
    diff(epsilont);
    diff(deltat);
    diff(dNt);
    diff(dTt)
];

disp('Completed reflectivity matrix construction.');

% Final synthetic gather
trueseis = G * real_R;
obseis_nref_ntr = reshape(trueseis, nref, ntr);

ntheta = 21;  % Reference angle
obseis_OA = zeros(nref, nazi);
for i = 1:nazi
    obseis_OA(:, i) = -obseis_nref_ntr(:, (i-1)*nangle + ntheta);
end

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

