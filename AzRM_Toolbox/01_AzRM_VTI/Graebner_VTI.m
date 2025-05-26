function [Rpp] = Graebner_VTI(angle, rho, C_VTI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact plane-wave reflection coefficients for VTI media
%
% Input:
%   angle : incidence angles (degrees)
%   rho   : density (g/cm^3)
%   C_VTI : stiffness matrix [6x6xN] for N layers
%
% Output:
%   Rpp   : PP reflection coefficients
%
% Reference:
%   Graebner, M. (1992). Geophysics, 57(11), 1512–1519.
%
%   Written by Rui Yang, Tongji University, 2024/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract stiffness terms
A = squeeze(C_VTI(1, 1, :));
C = squeeze(C_VTI(3, 3, :));
F = squeeze(C_VTI(1, 3, :));
L = squeeze(C_VTI(5, 5, :));

vp       = sqrt(C ./ rho);        % Vertical P-wave velocity
nLayers  = length(A);
nAngles  = length(angle);

Rpp      = zeros(nLayers - 1, nAngles);
pslowness = zeros(nLayers, nAngles);
sslowness = zeros(nLayers, nAngles);

for i = 1:nAngles
    theta = angle(i) * pi / 180;
    uC    = sin(theta) / vp(1);
    uuC   = uC^2;

    % ----------------------------------------------
    % Compute vertical slownesses (Graebner, Eq. 8-9)
    % ----------------------------------------------
    aux = A ./ L + L ./ C - ((F + L).^2) ./ (C .* L);
    K1  = rho ./ C + rho ./ L - aux * uuC;
    K2  = A ./ C * uuC - rho ./ C;
    K3  = uuC - rho ./ L;

    a = 1/sqrt(2) * sqrt(K1 - sqrt(K1.^2 - 4 * K2 .* K3));  % P-wave
    b = 1/sqrt(2) * sqrt(K1 + sqrt(K1.^2 - 4 * K2 .* K3));  % S-wave

    pslowness(:, i) = a;
    sslowness(:, i) = b;

    % ----------------------------------------------
    % Precompute auxiliary terms
    % ----------------------------------------------
    Ar = A ./ rho;
    Cr = C ./ rho;
    Lr = L ./ rho;

    aux1p = Ar * uuC + Lr .* a.^2 - 1;
    aux2p = Cr .* a.^2 + Lr * uuC - 1;
    aux1s = Ar * uuC + Lr .* b.^2 - 1;
    aux2s = Cr .* b.^2 + Lr * uuC - 1;

    for n = 1:nLayers - 1
        i1 = n;
        i2 = n + 1;

        % ------------------------------------------
        % Compute P and S eigenvector components
        % ------------------------------------------
        % Layer i1
        lalpha(i1, i) = sqrt(aux2p(i1) / (aux1p(i1) + aux2p(i1)));
        malpha(i1, i) = sqrt(aux1p(i1) / (aux1p(i1) + aux2p(i1)));

        lbeta(i1, i)  = sqrt(aux1s(i1) / (aux1s(i1) + aux2s(i1)));
        mbeta(i1, i)  = sqrt(aux2s(i1) / (aux1s(i1) + aux2s(i1)));

        % Layer i2
        lalpha(i2, i) = sqrt(aux2p(i2) / (aux1p(i2) + aux2p(i2)));
        malpha(i2, i) = sqrt(aux1p(i2) / (aux1p(i2) + aux2p(i2)));

        lbeta(i2, i)  = sqrt(aux1s(i2) / (aux1s(i2) + aux2s(i2)));
        mbeta(i2, i)  = sqrt(aux2s(i2) / (aux1s(i2) + aux2s(i2)));

        % ------------------------------------------
        % Build Q and P matrices (Graebner Eq. 13–14)
        % ------------------------------------------
        a1 = L(i1) * (a(i1) * lalpha(i1, i) + uC * malpha(i1, i));
        a2 = L(i2) * (a(i2) * lalpha(i2, i) + uC * malpha(i2, i));
        b1 = L(i1) * (b(i1) * mbeta(i1, i) - uC * lbeta(i1, i));
        b2 = L(i2) * (b(i2) * mbeta(i2, i) - uC * lbeta(i2, i));

        c1 = uC * lalpha(i1, i) * F(i1) + a(i1) * malpha(i1, i) * C(i1);
        c2 = uC * lalpha(i2, i) * F(i2) + a(i2) * malpha(i2, i) * C(i2);

        d1 = uC * mbeta(i1, i) * F(i1) - b(i1) * lbeta(i1, i) * C(i1);
        d2 = uC * mbeta(i2, i) * F(i2) - b(i2) * lbeta(i2, i) * C(i2);

        % Matrices P and Q
        P1 = [lalpha(i1, i), mbeta(i1, i); c1, d1];
        P2 = [lalpha(i2, i), mbeta(i2, i); c2, d2];

        Q1 = [malpha(i1, i), -lbeta(i1, i); a1, b1];
        Q2 = [malpha(i2, i), -lbeta(i2, i); a2, b2];

        % ------------------------------------------
        % Matrix inversion (manual or inv)
        % ------------------------------------------
        P1I = inv(P1);
        P2I = inv(P2);

        QP  = Q1 * P1I + Q2 * P2I;
        QPI = inv(QP);

        % ------------------------------------------
        % Reflection and transmission coefficients
        % ------------------------------------------
        Rd = P1I * QPI * (Q1 * P1I - Q2 * P2I) * P1;
        Td = 2 * P2I * QPI * Q1;

        Ru = P2I * QPI * (Q2 * P2I - Q1 * P1I) * P2;
        Tu = 2 * P1I * QPI * Q2;

        % Store PP reflection coefficient (real part only)
        Rpp(n, i) = real(Rd(1, 1));
    end
end
end
