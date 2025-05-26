%**************************************************************************
% Main Program for Azimuthally Anisotropic Reflectivity (AzRM) toolbox in VTI Media
%
%   Written by Rui Yang, Tongji University, 2024/06
%**************************************************************************

clc;
clear;
close all;

%% Load model and parameters
model      = load('loadmodel_VTI.txt');
parameters = load('parameters_VTI_model.txt');

thickness = model(:, 1);     % Thickness per layer
rho       = model(:, 4);     % Density

% Modeling parameters
Tmax     = parameters(1);    % Max modeling time (s)
delta_T  = parameters(2);    % Time interval (ms)
f1       = parameters(3);    % Start frequency (Hz)
f2       = parameters(4);    % End frequency (Hz)
theta1   = parameters(5);    % Start incidence angle (deg)
theta2   = parameters(6);    % End incidence angle (deg)
fdom     = parameters(9);    % Dominant frequency (Hz)

% Sampling and angle parameters
t        = (0:delta_T:Tmax)';
nSamples = round(Tmax / delta_T) + 1;
delta_F  = 1 / Tmax;
nF       = round((f2 - f1) / delta_F);

delta_theta = 1;
nangle      = round((theta2 - theta1) / delta_theta) + 1;
angle       = theta1:delta_theta:theta2;

% Frequency-domain Ricker wavelet settings
dt   = 0.002;    % Time interval of wavelet (s)
tlen = 0.3;      % Duration of wavelet (s)

%% Azimuthally Anisotropic Reflectivity Method (exact)
[uZ, c_layer] = AzRM_VTI(model, parameters);

% Display synthetic gather
figure('Name', 'Synthetic gathers in a VTI medium');
plotseis_AzRM(uZ, t, angle, [], 3);
ylabel('Time (s)');
xlabel('\theta (\circ)');
box on;
bigfont(gcf, 1.25, 1);

%% Ruger Approximation for VTI
obseis_Ruger = Ruger_VTI(Tmax, thickness, dt, fdom, tlen, angle, rho, c_layer);

figure('Name', 'Ruger method');
plotseis_AzRM(obseis_Ruger, t(1:end-1), angle, [], 3);
ylabel('Time (s)');
xlabel('\theta (\circ)');
box on;
bigfont(gcf, 1.25, 1);

%% Graebner Exact VTI Reflection Coefficient
Rpp_Graebner = Graebner_VTI(angle, rho, c_layer);

%% Automatic AVA Layer Detection
first_col = uZ(1:40, 1);
[~, sorted_first] = sort(first_col, 'descend');
first_layer = sorted_first(1);

second_start = first_layer + 10;
second_col   = uZ(second_start:end, 1);
[~, sorted_indices] = sort(abs(second_col), 'descend');
second_layer = second_start + sorted_indices(1) - 1;

%% AVA Analysis - First Layer
target = first_layer;
layer  = 1;

if uZ(target, 1) == 0
    error('The first element of uZ cannot be zero');
end

scaleFactor_RM = obseis_Ruger(target, 1) / uZ(target, 1);
uZ_new = uZ * scaleFactor_RM;

scaleFactor_Gr = obseis_Ruger(target, 1) / Rpp_Graebner(layer, 1);
Rpp_Graebner   = Rpp_Graebner * scaleFactor_Gr;

figure('Name', 'VTI Reflectivity AVO analysis - First Layer');
plot(angle, uZ_new(target, :), '-', 'Color', '#D95319', 'LineWidth', 3);
hold on;
plot(angle, obseis_Ruger(target, :), '-', 'Color', '#0072BD', 'LineWidth', 3);
plot(angle, Rpp_Graebner(layer, :), '--', 'Color', '#000000', 'LineWidth', 3);
text(-6.6, 0.152, 'a)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('P-wave reflection coefficient');
xlabel('\theta (\circ)');
legend('Reflectivity', 'Approx VTI', 'Exact VTI', 'Location', 'Northwest', 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
ylim([0.06 0.15]);

%% AVA Analysis - Second Layer
target = second_layer;
layer  = 2;

if uZ(first_layer, 1) == 0
    error('The first element of uZ cannot be zero');
end

scaleFactor_RM = obseis_Ruger(first_layer, 1) / uZ(first_layer, 1);
uZ_new = uZ * scaleFactor_RM;

scaleFactor_Gr = obseis_Ruger(target, 1) / Rpp_Graebner(layer, 1);
Rpp_Graebner   = Rpp_Graebner * scaleFactor_Gr;

figure('Name', 'VTI Reflectivity AVO analysis - Second Layer');
plot(angle, uZ_new(target, :), '-', 'Color', '#D95319', 'LineWidth', 3);
hold on;
plot(angle, obseis_Ruger(target, :), '-', 'Color', '#0072BD', 'LineWidth', 3);
plot(angle, Rpp_Graebner(layer, :), '--', 'Color', '#000000', 'LineWidth', 3);
text(-7.5, -0.057, 'b)', 'FontName', 'Times New Roman', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('P-wave reflection coefficient');
xlabel('\theta (\circ)');
legend('Reflectivity', 'Approx VTI', 'Exact VTI', 'Location', 'Northwest', 'FontSize', 16);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
ylim([-0.15 -0.06]);
