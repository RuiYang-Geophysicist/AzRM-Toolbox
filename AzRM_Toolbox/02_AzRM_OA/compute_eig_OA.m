function [D, T] = compute_eig_OA(px, py, rho, c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes eigenvalues and eigenvectors for orthorhombic media
% based on Fryer & Frazer (1987)
%
% INPUT:
%   px, py : horizontal slownesses
%   rho    : density
%   c      : stiffness matrix for a single layer (6x6)
%
% OUTPUT:
%   D      : Eigenvector matrix (6x6)
%   T      : Diagonal matrix of vertical slowness eigenvalues
%
% Reference:
%   Fryer & Frazer, GJI (1987), Eqns 3.9–3.11 and 4.2
%
%   Written by Rui Yang, Tongji University, 2024/06
%              Neng Lu, Jilin University, 2015/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------ Extract stiffness components ------------------
c11 = c(1,1); c12 = c(1,2); c13 = c(1,3);
c22 = c(2,2); c33 = c(3,3); c23 = c(2,3);
c44 = c(4,4); c55 = c(5,5); c66 = c(6,6);

b = zeros(6, 6);   % Eigenvector matrix
D = zeros(6, 6);
T = zeros(6, 6);

% ------------------------ Build system matrix A (Eq. 4.2) ----------------
s11 = px^2 * (c11 - c13^2 / c33) + py^2 * c66 - rho;
s12 = px * py * (c12 + c66 - c13 * c23 / c33);
s22 = px^2 * c66 + py^2 * (c22 - c23^2 / c33) - rho;

A1 = [0 0 -px*c13/c33 s11 s12 0]';  
A2 = [0 0 -py*c23/c33 s12 s22 0]';
A3 = [-px -py 0 0 0 -rho]';
A4 = [-1/c55 0 0 0 0 -px]';
A5 = [0 -1/c44 0 0 0 -py]';
A6 = [0 0 -1/c33 -px*c13/c33 -py*c23/c33 0]';
A  = [A1 A2 A3 A4 A5 A6];


% ------------------------ Solve for eigenvalues -------------------------
[~, eigvals] = eig(A);
qi = diag(eigvals);
sorted_qi = sort(qi, 'descend');

q6 = sorted_qi(1);  % qD_T
q5 = sorted_qi(2);  % qD_S
q4 = sorted_qi(3);  % qD_P
q3 = sorted_qi(6);  % qU_T
q2 = sorted_qi(5);  % qU_S
q1 = sorted_qi(4);  % qU_P

q = [q1; q2; q3; q4; q5; q6];

% ------------------------ Eigenvector construction ----------------------
for i = 1:6
    % Define auxiliary variables (Eqns 3.9–3.11)
    alpha = c11 * px^2 + c66 * py^2 + c55 * q(i)^2;
    beta  = c66 * px^2 + c22 * py^2 + c44 * q(i)^2;
    gamma = c55 * px^2 + c44 * py^2 + c33 * q(i)^2;
    delta = (c12 + c66) * px * py;
    eta   = (c13 + c55) * px * q(i);
    zeta  = (c23 + c44) * py * q(i);

    % ---------------- qP wave (Eq. 3.9) ----------------
    if i == 1 || i == 4
        ux = delta * zeta - eta * (beta - rho);
        uy = delta * eta - zeta * (alpha - rho);
        uz = (alpha - rho) * (beta - rho) - delta^2;

        tau1 = -c55 * (px * uz + q(i) * ux);
        tau2 = -c44 * (py * uz + q(i) * uy);
        tau3 = -c13 * px * ux - c23 * py * uy - c33 * q(i) * uz;

        uP   = [ux; uy; uz];
        tauP = [tau1; tau2; tau3];
        epsilon_p = 1/sqrt(abs(tauP'*uP+uP'*tauP));
        b(:,i) = epsilon_p.*[uP;tauP];

    % ---------------- qS1 wave (Eq. 3.10) ----------------
    elseif i == 2 || i == 5
        ux = eta * zeta - delta * (gamma - rho);
        uy = (alpha - rho) * (gamma - rho) - eta^2;
        uz = eta * delta - zeta * (alpha - rho);

        tau1 = -c55 * (px * uz + q(i) * ux);
        tau2 = -c44 * (py * uz + q(i) * uy);
        tau3 = -c13 * px * ux - c23 * py * uy - c33 * q(i) * uz;

        uSV   = [ux; uy; uz];
        tauSV = [tau1; tau2; tau3];
        epsilon_SV = 1/sqrt(abs(tauSV'*uSV+uSV'*tauSV));
        b(:,i) = epsilon_SV.*[uSV;tauSV];

        % Fallback if uSV is zero
        if norm(uSV) == 0
            ux = (beta - rho) * (gamma - rho) - zeta^2;
            uy = eta * zeta - delta * (gamma - rho);
            uz = delta * zeta - eta * (beta - rho);

            tau1 = -c55 * (px * uz + q(i) * ux);
            tau2 = -c44 * (py * uz + q(i) * uy);
            tau3 = -c13 * px * ux - c23 * py * uy - c33 * q(i) * uz;

            uSV   = [ux; uy; uz];
            tauSV = [tau1; tau2; tau3];
            epsilon_SV = 1/sqrt(abs(tauSV'*uSV+uSV'*tauSV));
            b(:,i) = epsilon_SV.*[uSV;tauSV];
        end

    % ---------------- qS2 wave (Eq. 3.11) ----------------
    elseif i == 3 || i == 6
        ux = (beta - rho) * (gamma - rho) - zeta^2;
        uy = eta * zeta - delta * (gamma - rho);
        uz = delta * zeta - eta * (beta - rho);

        tau1 = -c55 * (px * uz + q(i) * ux);
        tau2 = -c44 * (py * uz + q(i) * uy);
        tau3 = -c13 * px * ux - c23 * py * uy - c33 * q(i) * uz;

        uSH   = [ux; uy; uz];
        tauSH = [tau1; tau2; tau3];
        epsilon_SH = 1/sqrt(abs(tauSH'*uSH+uSH'*tauSH));
        b(:,i) = epsilon_SH.*[uSH;tauSH];

        % Fallback if uSH is zero
        if norm(uSH) == 0
            ux = eta * zeta - delta * (gamma - rho);
            uy = (alpha - rho) * (gamma - rho) - eta^2;
            uz = eta * delta - zeta * (alpha - rho);

            tau1 = -c55 * (px * uz + q(i) * ux);
            tau2 = -c44 * (py * uz + q(i) * uy);
            tau3 = -c13 * px * ux - c23 * py * uy - c33 * q(i) * uz;

            uSH   = [ux; uy; uz];
            tauSH = [tau1; tau2; tau3];
            epsilon_SH = 1/sqrt(abs(tauSH'*uSH+uSH'*tauSH));
            b(:,i) = epsilon_SH.*[uSH;tauSH];
        end
    end
end

% Output assignments
D = b;
T = diag(q);

end
