function [uZ] = AzRM_OA(model, parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuthally Anisotropic Reflectivity Method in Orthorhombic Media
%
% Elastic Model:
%   Each row of `model` must be formatted as:
%     [thickness (km), vp (km/s), vs (km/s), rho (g/cm^3), epsilon, delta, gamma, e]
%
% Parameters:
%   Inputs:
%     Tmax     : Maximum modeling time (s)
%     delta_T  : Time sampling interval (ms)
%     f1, f2   : Frequency range (Hz)
%     theta1, theta2 : Incidence angle range (deg)
%     phi1, phi2     : Azimuth angle range (deg)
%     fdom     : Dominant frequency of Ricker wavelet (Hz)
%     phi0     : Fracture orientation (deg)
%
% Output:
%   uZ : Synthetic PP-wave gathers in time domain [samples × angles × azimuths]
%
% References:
%   Fryer & Frazer (1984, 1987), Mallick & Frazer (1988)
%
% Note:
%   This implementation assumes idealized VTI-to-orthorhombic symmetry 
%   transformation. Use with caution.
%
%   Written by Rui Yang, Tongji University, 2024/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------- Model Extraction ----------------------------
thickness = model(:, 1);         % Layer thicknesses
rho       = model(:, 4);         % Density vector
c_layer   = Stiffness_matrix_OA(model);  % Orthorhombic stiffness tensors

% ------------------------ Parameter Extraction -------------------------
Tmax     = parameters(1);
delta_T  = parameters(2);
f1       = parameters(3);
f2       = parameters(4);
theta1   = parameters(5);
theta2   = parameters(6);
phi1     = parameters(7);
phi2     = parameters(8);
fdom     = parameters(9);
phi0     = parameters(10);

% ------------------------ Time/Frequency/Azimuth Grid -------------------
t        = (0:delta_T:Tmax)';
nSamples = round(Tmax / delta_T) + 1;
delta_F  = 1 / Tmax;
nF       = round((f2 - f1) / delta_F);

delta_theta = 1;
angle       = theta1:delta_theta:theta2;
nangle      = length(angle);

delta_phi = 10;
nazimuth = round((phi2 - phi1) / delta_phi)+1;        % # of azimuth to model 

% ------------------------ Preallocate Outputs ---------------------------
Rw_pp   = zeros(nF, nangle, nazimuth);
Wavelet=zeros(1,1);
uw      = zeros(nF, nangle, nazimuth);
uZ      = zeros(nSamples, nangle, nazimuth);

% ------------------------ Frequency Loop -------------------------------
iF = 1;
nwrite = 10;
tstart = clock;

for f = f1:delta_F:f2
    if(rem(iF,nwrite)==0)
        tnow = clock;
        timeused = etime(tnow, tstart);
        timeperfreq = timeused / iF;
        timeremaining = timeperfreq * (nF - iF);
        disp(['Computation at frequency ', num2str(f), ' Hz']);
        disp(['Time used: ', int2str(timeused), 's | Time remaining: ', int2str(timeremaining), 's']);
    end

    iphi = 1;
    for phi=phi1:delta_phi:phi2
        azimuth_shifted = phi + phi0;

        [Rpp, ~, ~] = compute_Rpp_OA(f, angle, thickness, c_layer, rho, azimuth_shifted);

        Rw_pp(iF, :, iphi) = Rpp;
        iphi = iphi + 1;
    end

    % Compute Ricker wavelet amplitude (frequency domain)
    Wavelet(iF,1) = 2 / sqrt(pi) .* (f / fdom).^2 * exp(-(f / fdom).^2);

    % Modulate reflectivity with wavelet
    uw(iF, :, :) = Rw_pp(iF, :, :) * Wavelet(iF);

    iF = iF + 1;
end

% ------------------------ Inverse FFT to Time Domain -------------------
for azi = 1:nazimuth
    for ang = 1:nangle
        uZ(:, ang, azi) = real(ifft(uw(:, ang, azi), nSamples));
    end
end

disp('Completed Azimuthal Anisotropic Reflectivity');

end
