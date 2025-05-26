function [C_OA] = Stiffness_matrix_OA(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute stiffness tensor of orthorhombic medium with vertical fractures
% using effective medium theory (VTI base + fracture perturbation)
%
% INPUT:
%   model : [n × 8] array with columns:
%       [1] thickness (km)
%       [2] vp (km/s)
%       [3] vs (km/s)
%       [4] rho (g/cm³)
%       [5] epsilon
%       [6] delta
%       [7] gamma
%       [8] vertical fracture density (e)
%
% OUTPUT:
%   C_OA  : [6 × 6 × nL] orthorhombic stiffness tensor (in GPa)
%
% REFERENCES:
%   Thomsen (1986), Geophysics, VTI anisotropy
%   Schoenberg & Douma (1988), T-matrix fracture weaknesses
%   Schoenberg & Helbig (1997), Orthorhombic media from vertical fractures
%
%   Written by Rui Yang, Tongji University, 2024/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------- Extract model parameters -----------------
vp      = model(:,2);             % P-wave velocity (km/s)
vs      = model(:,3);             % S-wave velocity (km/s)
rho     = model(:,4);             % Density (g/cm^3)
epsilon = model(:,5);             % Thomsen parameter
delta   = model(:,6);             % Thomsen parameter
gamma   = model(:,7);             % Thomsen parameter
e       = model(:,8);             % Fracture density
nL      = size(model,1);          % Number of layers

% ----------------- Compute VTI baseline stiffness -----------------
c33b = (vp.^2) .* rho;                        % C33 baseline
c55b = (vs.^2) .* rho;                        % C55 baseline
c44b = c55b;                                  % C44 baseline
c11b = c33b .* (1 + 2 .* epsilon);            % C11 from epsilon
c13b = sqrt(2 .* delta .* c33b .* (c33b - c55b) + (c33b - c55b).^2) - c55b;
c12b = c11b - 2 .* c55b;                      % C12 from C11, C55
c66b = c44b .* (1 + 2 .* gamma);              % C66 from gamma

% ----------------- Compute fracture weaknesses -----------------
g  = (vs ./ vp).^2;
dN = 4 .* e ./ (3 .* g .* (1 - g));           % Normal weakness
dT = 16 .* e ./ (3 .* (3 - 2 .* g));          % Tangential weakness

% ----------------- Compute fracture perturbations -----------------
c11f = -c11b .* dN;
c12f = -c12b .* dN;
c13f = -c13b .* dN;
c22f = -c12b.^2 ./ c11b .* dN;
c23f = -c12b .* c13b ./ c11b .* dN;
c33f = -c13b.^2 ./ c11b .* dN;
c55f = -c55b .* dT;
c66f = -c66b .* dT;

% ----------------- Build final stiffness tensors -----------------
C_OA  = zeros(6,6,nL);  % Final orthorhombic stiffness tensor
C_VTI = zeros(6,6,nL);  % VTI baseline
C_f   = zeros(6,6,nL);  % Fracture perturbation

for c = 1:nL
    % VTI matrix
    C_VTI(:,:,c) = [c11b(c) c12b(c) c13b(c) 0        0        0;
                    c12b(c) c11b(c) c13b(c) 0        0        0;
                    c13b(c) c13b(c) c33b(c) 0        0        0;
                    0       0       0       c55b(c) 0        0;
                    0       0       0       0       c55b(c) 0;
                    0       0       0       0       0       c66b(c)];

    % Fracture correction matrix
    C_f(:,:,c) = [c11f(c) c12f(c) c13f(c) 0        0        0;
                  c12f(c) c22f(c) c23f(c) 0        0        0;
                  c13f(c) c23f(c) c33f(c) 0        0        0;
                  0       0       0       0        0        0;
                  0       0       0       0       c55f(c) 0;
                  0       0       0       0        0       c66f(c)];

    % Total stiffness = VTI + fracture perturbation
    C_OA(:,:,c) = C_VTI(:,:,c) + C_f(:,:,c);
end

end
