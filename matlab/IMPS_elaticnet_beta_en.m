clear; clc; close all

% Load data from a text file
filename = '0.4 V.txt'; % Insert the file name
data = load(filename);

frequency = data(:, 1);
H_re = data(:, 2);
H_im = data(:, 3);
% frequency = flipud(frequency);

% Uncomment the following lines for interpolation if needed
% f = logspace(log10(frequency(1)), log10(frequency(end)), 1.2 * length(frequency));
% H_re_new = interp1(frequency, H_re, f);
% H_re_new = H_re_new';
% H_im_new = interp1(frequency, H_im, f);
% H_im_new = H_im_new';
% data_iut_new = [f' H_re_new H_im_new];
% data_iut_new = flipud(data_iut_new);

% H_re = flipud(H_re);
% H_im = flipud(H_im);
omega = 2 * pi * frequency;
tau = 1 ./ omega;

% Create a log-spaced tau array
tau = logspace(log10(tau(1)), log10(tau(end)), 10 * length(tau));

N = length(frequency);
M = length(tau);

% Initialize matrices A_re and A_im
A_re = zeros(N, M);
A_im = zeros(N, M);
beta = 1;
for j = 1:N
    for k = 1:M
        A_re(j, k) = real(1./(1 + 1i * omega(j) * tau(k)).^beta);
        A_im(j, k) = imag(1./(1 + 1i * omega(j) * tau(k)).^beta);
        % data_cell = 1./(1+1i*omega*tau_cell).^alpha(j);
        % cell_re = real(data_cell);
        % cell_im = imag(data_cell);
    end
end

% Define Elastic Net hyperparameters
% When lambda is small, regularization is weak, and the model fits the training data more easily but may lead to overfitting.
% When lambda is large, regularization is strong, and the model tends to be sparse (feature selection), reducing the risk of overfitting.
lambda = 0.5;
alpha = 0.5;  % You can adjust this parameter between 0 and 1 for different mixing ratios.

% Define the optimization problem using CVX with Elastic Net regularization
cvx_begin
    variable x(M)
    minimize(norm(H_re - A_re * x, 2) + norm(H_im - A_im * x, 2) + ...
        lambda * ((1 - alpha) * norm(x, 2) / length(H_re) + alpha * norm(x, 1) / length(H_re)))
cvx_end

% The optimal solution is stored in the variable 'x'
g = x;

%% Reconstruct real and imaginary parts using g(tau)
H_re_fit = zeros(N, 1);
H_im_fit = zeros(N, 1);

for j = 1:N
    for k = 1:length(tau)
        H_re_fit(j) = H_re_fit(j) + real((g(k) / (1 + 1j * omega(j) * tau(k))));
        H_im_fit(j) = H_im_fit(j) + imag((g(k) / (1 + 1j * omega(j) * tau(k))));
    end
end

%% Calculate Gartner current, recombination current, and total photocurrent
g_e = 0;
g_Gar = 0;
g_ph = 0;
g_w = 0;

for j = 1:length(tau)
    if g(j) > 0
        g_e = g_e + g(j); % recombination current
        g_w = g_w + g(j) * tau(j); % weighted average for calculating kt_mean and kr_mean
    else
        g_Gar = g_Gar + g(j); % Gartner current
    end
end

g_ph = g_Gar + g_e; % photocurrent
tau_mean = g_w / g_e;
kt_mean = (g_ph / tau_mean) / g_Gar;
kr_mean = (1 / tau_mean) - kt_mean;

% Save data
figure(2)
subplot(2, 2, [1, 3]);
plot(H_re, H_im, '-o', H_re_fit, H_im_fit, '-*', 'linewidth', 2)
xlabel('Re(H)')
ylabel('Im(H)')
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha');
set(gca, 'FontSize', 16, 'Fontname', 'Times', 'linewidth', 1.5);

subplot(2, 2, 2);
semilogx(frequency, H_re, '-o')
hold on
semilogx(frequency, H_re_fit, '-*', 'linewidth', 2)
xlabel('Frequency (Hz)')
ylabel('Re(H)')
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha');
set(gca, 'FontSize', 16, 'Fontname', 'Times', 'linewidth', 1.5);

subplot(2, 2, 4);
semilogx(frequency, H_im, '-o')
hold on
semilogx(frequency, H_im_fit, '-*', 'linewidth', 2)
xlabel('Frequency (Hz)')
ylabel('Im(H)')
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha');
set(gca, 'FontSize', 16, 'Fontname', 'Times', 'linewidth', 1.5);

figure(3)
subplot(2, 1, 1)
semilogx(tau, g, 'linewidth', 2)
xlabel('\tau (s)')
ylabel('g(\tau)')
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha');
set(gca, 'FontSize', 16, 'Fontname', 'Times', 'linewidth', 1.5);
set(gca, 'xaxislocation', 'top', 'yaxislocation', 'left');

subplot(2, 1, 2)
freq_fit = 1./(2 * pi * tau);
freq_fit = flipud(freq_fit);
semilogx(freq_fit, g, 'linewidth', 2)
xlabel('Frequency (Hz)')
set(gca, 'XDir', 'reverse') % Reverse X direction
ylabel('g(\tau)')
grid on
set(gca, 'GridLineStyle', ':');
set(gca, 'GridAlpha');
set(gca, 'FontSize', 16, 'Fontname', 'Times', 'linewidth', 1.5);
set(gca, 'xaxislocation', 'bottom', 'yaxislocation', 'left');

%% Save data
Data_H = [frequency, H_re, H_im, H_re_fit, H_im_fit];
file1 = '_H';
newFilename1 = [filename, file1, '.txt'];
dlmwrite(newFilename1, Data_H, 'delimiter', '\t');

Data_DRT = [freq_fit', tau', g];
file2 = '_DRT';
newFilename2 = [filename, file2, '.txt'];
dlmwrite(newFilename2, Data_DRT, 'delimiter', '\t');

yita = kt_mean / (kr_mean + kt_mean);
Data_JV = [g_e, g_Gar, g_ph, kr_mean, kt_mean, yita];
file3 = '_JV';
newFilename3 = [filename, file3, '.txt'];
dlmwrite(newFilename3, Data_JV, 'delimiter', '\t');