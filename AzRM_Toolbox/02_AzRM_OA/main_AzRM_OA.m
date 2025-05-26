%**************************************************************************
% Main Program for Azimuthally Anisotropic Reflectivity in Orthorhombic Media
%
%   Written by Rui Yang, Tongji University, 2024/06
%**************************************************************************

clc;
clear;
close all;

%% -------------------------- Load Model and Parameters --------------------------
model      = load('loadmodel_OA.txt');         % Input elastic model
parameters = load('parameters_OA_model.txt');  % Simulation parameters
thickness  = model(:, 1);                      % Layer thicknesses (km)

% Extract parameters
Tmax    = parameters(1);    %maximum modeling time (sec);
delta_T = parameters(2);    %time interval (msec);
f1      = parameters(3);    %initial frequency (Hz) to model;
f2      = parameters(4);    %final frequency (Hz) to model;
theta1  = parameters(5);    %initial incidence to model (degree);
theta2  = parameters(6);    %final incidence to model (degree);
phi1    = parameters(7);    %initial azimuth to model (degree);
phi2    = parameters(8);    %final azimuth to model (degree);
fdom    = parameters(9);    %the dominant frequency of the wavelet;
phi0    = parameters(10);   %fracture orientation;

% Sampling grid
t        = (0:delta_T:Tmax)';
nSamples = round(Tmax / delta_T) + 1;
delta_F  = 1 / Tmax;
nF       = round((f2 - f1) / delta_F);

% Incidence and azimuth angles
delta_theta = 1;
angle       = theta1:delta_theta:theta2;
nangle      = length(angle);

delta_phi = 10;
azimuth   = phi1:delta_phi:phi2;
nazimuth  = length(azimuth);

% Wavelet parameters
dt   = 0.002;  % s
tlen = 0.3;    % s

%% -------------------------- AzRM Modeling --------------------------
disp('Start AzRM method');
uZ = AzRM_OA(model, parameters);

%% -------------------------- Ruger Approximation --------------------------
disp('Start Ruger method (need some time because of depth-time conversion)');
tic;
obseis_OA = Ruger_OA(Tmax, thickness, dt, fdom, tlen, angle, model, phi1, phi2, delta_phi, phi0);
t_elapsed = toc;
disp(['Time used: ', num2str(t_elapsed), ' s']);
disp('Completed Ruger method');

% Display Ruger gather
figure('Name','Ruger method');
plotseis_AzRM(obseis_OA, t(1:end-1), azimuth, [], 2);
ylabel('Time (s)'); xlabel('\phi (\circ)');
set(gca, 'XTick', 0:60:360);
box on;
bigfont(gcf, 1.25, 1);

%% -------------------------- AVA Analysis at Fixed Azimuth --------------------------
nphi = 1;  % Select azimuth index

% Display angle gather at selected azimuth
figure('Name','Angle gathers at fixed azimuth');
plotseis_AzRM(uZ(:,:,nphi), t, angle, [], 2);
ylabel('Time (s)'); xlabel('\theta (\circ)');
box on;
bigfont(gcf, 1.25, 1);

% Identify top and bottom reflectors (based on amplitude)
first_col = uZ(1:40, 1, nphi);
[~, sorted_first] = sort(abs(first_col), 'descend');
first_layer1 = sorted_first(1);

second_start = first_layer1 + 10; %Notice: it is possible that there are larger values than second layer around first layer
second_col = uZ(second_start:end, 1, nphi);
[~, sorted_indices] = sort(abs(second_col), 'descend');
second_layer1 = second_start + sorted_indices(1) - 1;

% Plot AVA curves
figure('Name','AVA at fixed azimuth');
plot(angle, 100 * uZ(first_layer1,:,nphi), 'Color','#0072BD', 'LineWidth', 3); hold on;
plot(angle, 100 * uZ(second_layer1,:,nphi), 'Color','#D95319', 'LineWidth', 3);
text(-4.5, 0.3, 'a)', 'FontName','Times New Roman','FontSize', 24, 'FontWeight','bold');
ylabel('Amplitude'); xlabel('\theta (\circ)');
legend('Top layer', 'Bottom layer', 'Location', 'Best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
ylim([-0.3, 0.3]);

%% -------------------------- AVAZ Analysis at Fixed Angle --------------------------
ntheta = 31;  % Select incidence angle index
uZ_AVAZ = reshape(uZ(:, ntheta, :), [nSamples, nazimuth]);

% Display azimuthal gather
figure('Name','Azimuthal gathers at fixed angle');
plotseis_AzRM(uZ_AVAZ, t, azimuth, [], 2);
ylabel('Time (s)'); xlabel('\phi (\circ)');
set(gca, 'XTick', 0:60:360);
box on;
bigfont(gcf, 1.25, 1);

% Identify reflection layers
first_col = uZ_AVAZ(1:40, 1);
[~, sorted_first] = sort(abs(first_col), 'descend');
first_layer2 = sorted_first(1);

second_start = first_layer2 + 10;
second_col = uZ_AVAZ(second_start:end, 1);
[~, sorted_indices] = sort(second_col, 'ascend');
second_layer2 = second_start + sorted_indices(1) - 1;

% Plot AVAZ curves (AzRM)
figure('Name','AVAZ Reflectivity (AzRM)');
plot(azimuth, 100 * uZ_AVAZ(first_layer2, :), '-', 'Color','#0072BD', 'LineWidth', 2.5); hold on;
plot(azimuth, 100 * uZ_AVAZ(second_layer2, :), '-', 'Color','#D95319', 'LineWidth', 2.5); 
text(-67, 0.005, 'a)', 'FontName','Times New Roman','FontSize', 24, 'FontWeight','bold');
ylabel('P-wave reflection coefficient'); xlabel('\phi (\circ)');
legend('Top layer', 'Bottom layer', 'Location', 'Best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'XTick', 0:60:360);
ylim([-0.13 0]);

% Plot AVAZ curves (Ruger)
figure('Name','AVAZ Reflectivity (Ruger)');
plot(azimuth, obseis_OA(first_layer2 - 1, :), '-', 'Color','#0072BD', 'LineWidth', 2.5); hold on;
plot(azimuth, obseis_OA(second_layer2 - 9, :), '-', 'Color','#D95319', 'LineWidth', 2.5);
text(-68, 0.045, 'b)', 'FontName','Times New Roman','FontSize', 24, 'FontWeight','bold');
ylabel('P-wave reflection coefficient'); xlabel('\phi (\circ)');
legend('Top layer', 'Bottom layer', 'Location', 'Best');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'XTick', 0:60:360);
ylim([-0.03 0.042]);

%% -------------------------- Polar Plot for Fracture Orientation --------------------------
amplitude2    = 1000 * uZ_AVAZ(second_layer2, :);
azimuth_rad   = deg2rad([azimuth, azimuth(1)]);
amplitude2    = [amplitude2, amplitude2(1)];

figure('Name','Fracture Orientation Polar Plot');
polarplot(azimuth_rad, amplitude2, 'Color','#0072BD', 'LineWidth', 3);
text(2.5, 0.75, 'd)', 'FontName','Times New Roman','FontSize', 24, 'FontWeight','bold');
ax = gca;
ax.ThetaTick = 0:30:330;
ax.ThetaTickLabel = arrayfun(@(x) sprintf('%d°', x), ax.ThetaTick, 'UniformOutput', false);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
