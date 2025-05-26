function [Rpp, Rpsv, Rpsh] = compute_Rpp_VTI(f, angle, thickness, c_layer, den, phi)
% compute_Rpp_VTI computes PP, PSv, PSh reflection coefficients using
% Fryer and Frazer's eigenstructure-based method in VTI media.
%
%   Written by Rui Yang, Tongji University, 2024/06
%              Neng Lu, Jilin University, 2015/05

nL   = length(den);              % Number of layers
Na   = length(angle);            % Number of incidence angles
phi  = phi * pi / 180;           % Convert to radians

% Horizontal slowness (based on top layer vertical P velocity)
Sp = 1 / sqrt(c_layer(3, 3, 1) / den(1));

% Preallocate matrices
T0(:,:,:) = zeros(6,6,nL);  % initialize eigenvalue T0 and eigenvector D0 matrices 
D0(:,:,:) = zeros(6,6,nL);  % as waves propagating vertically
T(:,:,:)  = zeros(6,6,nL);  % initialize eigenvalue T and eigenvector D matrices
D(:,:,:)  = zeros(6,6,nL);  % as waves propagating arbitrarily
Ed(:,:,:) = zeros(3,3,nL-1);% initialize the phase shift matrix Ed,Eu
Eu(:,:,:) = zeros(3,3,nL-1);
ru(:,:,:) = zeros(3,3,nL);  % initialize the reflection matrices ru,rd of i layer
rd(:,:,:) = zeros(3,3,nL);
tu(:,:,:) = zeros(3,3,nL);  % initialize the transmission matrices tu,td of i layer
td(:,:,:) = zeros(3,3,nL);
Ru(:,:,:) = zeros(3,3,nL);  % initialize the reflection matrices Ru,Rd of i+1 layer
Rd(:,:,:) = zeros(3,3,nL);  
Tu(:,:,1) = eye(3);                % initialize the transmission matrices Tu,Rd of i+1 layer
Td(:,:,1) = eye(3);    

Q    = zeros(6, 6, nL);
Q11  = zeros(3, 3, nL); Q12 = zeros(3, 3, nL);
Q21  = zeros(3, 3, nL); Q22 = zeros(3, 3, nL);

Rpp  = zeros(1, Na);
Rpsv = zeros(1, Na);
Rpsh = zeros(1, Na);

% Angular frequency
w = 2 * pi * f;

%% Compute T0 (eigenvalues at normal incidence)
px = 0;
py = 0;
for ii = 1:nL
    [~, T0(:, :, ii)] = compute_eig_VTI(px, py, den(ii), c_layer(:, :, ii));
end

%% Loop over incidence angles
for i = 1:Na
    theta = angle(i) * pi / 180;
    px = Sp * sin(theta) * cos(phi);
    py = Sp * sin(theta) * sin(phi);

    % Compute eigenvalues/vectors for each layer
    for ii = 1:nL
        [D(:, :, ii), T(:, :, ii)] = compute_eig_VTI(px, py, den(ii), c_layer(:, :, ii));
    end

    %% Compute phase shift matrices (Eu, Ed)
    for ii = 1:nL-1
        h = thickness(ii);
        e = zeros(1, 6);
        for jj = 1:6
            e(jj) = exp(T0(jj, jj, ii) * h * w * 1i);
        end
        Eu(:, :, ii) = diag([e(1), e(2), e(3)]);
        Ed(:, :, ii) = diag([e(4), e(5), e(6)]);
    end

    %% Compute reflection/transmission at each interface
    for ii = 2:nL
        Q(:, :, ii)   = D(:, :, ii) \ D(:, :, ii-1);  % Fryer & Frazer (1984), Eq. (4.22)
        Q11(:, :, ii) = Q(1:3, 1:3, ii);
        Q12(:, :, ii) = Q(1:3, 4:6, ii);
        Q21(:, :, ii) = Q(4:6, 1:3, ii);
        Q22(:, :, ii) = Q(4:6, 4:6, ii);

        ru(:, :, ii)  = Q21(:, :, ii) / Q11(:, :, ii);
        rd(:, :, ii)  = -Q11(:, :, ii) \ Q12(:, :, ii);
        tu(:, :, ii)  = inv(Q11(:, :, ii));
        td(:, :, ii)  = Q22(:, :, ii) - Q21(:, :, ii) / Q11(:, :, ii) * Q12(:, :, ii);  % Eq. (4.7)
    end

    %% Recursive computation from free surface to last layer
    I = eye(3);
    for ii = 2:nL
        jj = ii - 1;

        Tu(:, :, ii) = Tu(:, :, jj) / (I - Eu(:, :, jj) \ rd(:, :, ii) * Ed(:, :, jj) * Ru(:, :, jj)) ...
                       / Eu(:, :, jj) * tu(:, :, ii);

        Rd(:, :, ii) = Rd(:, :, jj) + Tu(:, :, jj) / Eu(:, :, jj) * rd(:, :, ii) * Ed(:, :, jj) ...
                       / (I - Ru(:, :, jj) / Eu(:, :, jj) * rd(:, :, ii) * Ed(:, :, jj)) * Td(:, :, jj);

        Ru(:, :, ii) = ru(:, :, ii) + td(:, :, ii) * Ed(:, :, jj) * Ru(:, :, jj) ...
                       / (I - Eu(:, :, jj) \ rd(:, :, ii) * Ed(:, :, jj) * Ru(:, :, jj)) ...
                       / Eu(:, :, jj) * tu(:, :, ii);

        Td(:, :, ii) = td(:, :, ii) * Ed(:, :, jj) ...
                       / (I - Ru(:, :, jj) / Eu(:, :, jj) * rd(:, :, ii) * Ed(:, :, jj)) * Td(:, :, jj);
    end

    %% Extract reflection coefficients at surface
    Rpp(i)  = Rd(1, 1, nL)';
    Rpsv(i) = Rd(1, 2, nL)';
    Rpsh(i) = Rd(1, 3, nL)';
end

end
