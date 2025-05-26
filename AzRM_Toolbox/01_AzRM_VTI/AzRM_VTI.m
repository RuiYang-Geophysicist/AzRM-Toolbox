function [uZ, c_layer] = AzRM_VTI(model, parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuthally Anisotropic Reflectivity Method in VTI Media
%
% Elastic Model:
%   model: each row defines a layer:
%     [thickness(km), vp(km/s), vs(km/s), density(g/cm^3), epsilon, delta, gamma]
%
% Parameters:
%   Inputs:
%     T_max    : maximum modeling time (s)
%     delta_T  : time sampling interval (ms)
%     f1, f2   : initial/final frequency (Hz)
%     theta1, theta2 : incidence angle range (deg)
%     fdom     : dominant frequency of the Ricker wavelet (Hz)
%
% Outputs:
%     uZ       : PP-wave angle gathers in time domain
%     c_layer  : stiffness tensors of each layer
%
% References:
%   Fryer & Frazer, GJI 1984, 1987
%   Mallick & Frazer, Geophysics 1988
%
% NOTE:
%   This code is provided as-is without warranty. Use at your own risk.
%
%   Written by Rui Yang, Tongji University, 2024/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input model parameters
thickness = model(:, 1);          % Layer thicknesses
rho       = model(:, 4);          % Density

% Compute the stiffness matrix for VTI media
c_layer = Stiffness_matrix_VTI(model);

% Input modeling parameters
Tmax     = parameters(1);         % Max modeling time (s)
delta_T  = parameters(2);         % Time step (ms)
f1       = parameters(3);         % Start frequency (Hz)
f2       = parameters(4);         % End frequency (Hz)
theta1   = parameters(5);         % Min incidence angle (deg)
theta2   = parameters(6);         % Max incidence angle (deg)
fdom     = parameters(9);         % Dominant frequency (Hz)

% Set up time and frequency
t         = (0:delta_T:Tmax)';
nSamples  = round(Tmax / delta_T) + 1;
delta_F   = 1 / Tmax;
nF        = round((f2 - f1) / delta_F);

% Set up angle range
delta_theta = 1;
nangle      = round((theta2 - theta1) / delta_theta) + 1;
angle       = theta1:delta_theta:theta2;

% Preallocate output matrices
Rw_pp   = zeros(nF, nangle);
Wavelet = zeros(nF, 1);
uw      = zeros(nF, nangle);
uZ      = zeros(nSamples, nangle);

% Main loop over frequency
iF     = 1;
nwrite = 10;
tstart = clock;

for f = f1:delta_F:f2
    if(rem(iF,nwrite)==0)
        disp(['Computation started at frequency ' num2str(f)]);
        tnow          = clock;
        timeused      = etime(tnow, tstart);
        timeperfreq   = timeused / iF;
        timeremaining = timeperfreq * (nF - iF);
        disp(['Time used: ' int2str(timeused) 's, Remaining: ' int2str(timeremaining) 's']);
    end

    phi = 0;
    [Rpp, ~, ~] = compute_Rpp_VTI(f, angle, thickness, c_layer, rho, phi);

    % Compute the overall PP-wave reflection coefficient matrix of all layers for each frequency
    Rw_pp(iF, :) = Rpp;

    % Ricker wavelet in frequency domain
    Wavelet(iF) = 2 / sqrt(pi) * (f / fdom).^2 .* exp(-(f / fdom).^2);

    % Frequency-domain displacement
    uw(iF, :) = Rw_pp(iF, :) .* Wavelet(iF);

    iF = iF + 1;
end

% Convert to time domain using IFFT
for ang = 1:nangle
    uZ(:, ang) = real(ifft(uw(:, ang), nSamples));
end

% Final report
disp('Completed Azimuthal Anisotropic Reflectivity');
tnow     = clock;
timeused = etime(tnow, tstart);
disp(['Total time: ' int2str(timeused) 's']);

end
