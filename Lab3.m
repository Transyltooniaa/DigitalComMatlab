% Parameters
N = 1000; % Number of symbols
SNR_dB = 20; % Signal-to-Noise Ratio in dB
M = 4; % Modulation order (QPSK)
EbNo = 10^(SNR_dB/10); % Convert SNR from dB to linear scale

% Transmit data (random bits)
data = randi([0 M-1], 1, N); 

% QPSK Modulation
modulated_signal = pskmod(data, M, pi/M); 

% Rayleigh fading channel (complex Gaussian)
h = (randn(1, N) + 1i * randn(1, N)) / sqrt(2); 
received_signal = h .* modulated_signal; 

% Additive White Gaussian Noise (AWGN)
noise = (randn(1, N) + 1i * randn(1, N)) / sqrt(2*EbNo); 
received_signal_with_noise = received_signal + noise;

% Plot transmitted signal (real part)
figure;
subplot(3,1,1);
plot(real(modulated_signal), 'b');
title('Transmitted Signal (Real Part)');
xlabel('Symbol Index');
ylabel('Amplitude');
grid on;

% Plot fading channel (real and imaginary parts)
subplot(3,1,2);
plot(real(h), 'r');
hold on;
plot(imag(h), 'g');
title('Rayleigh Fading Channel (Real and Imaginary Parts)');
xlabel('Symbol Index');
ylabel('Amplitude');
legend('Real Part', 'Imaginary Part');
grid on;

% Plot received signal (real part after noise addition)
subplot(3,1,3);
plot(real(received_signal_with_noise), 'k');
title('Received Signal with Noise (Real Part)');
xlabel('Symbol Index');
ylabel('Amplitude');
grid on;

% Magnitude of received signal (Rayleigh distribution)
magnitude_received_signal = abs(received_signal);

% Plot histogram of the magnitude (empirical Rayleigh distribution)
figure;
histogram(magnitude_received_signal, 'Normalization', 'pdf');
hold on;

% Overlay theoretical Rayleigh distribution
x = linspace(0, max(magnitude_received_signal), 100);
sigma = 1/sqrt(2); % Rayleigh parameter
rayleigh_pdf = (x ./ (sigma^2)) .* exp(-(x.^2) / (2*sigma^2));
plot(x, rayleigh_pdf, 'r', 'LineWidth', 2);

title('Histogram of Received Signal Magnitude (Rayleigh Distribution)');
xlabel('Magnitude');
ylabel('Probability Density');
legend('Empirical Distribution', 'Theoretical Rayleigh PDF');
grid on;
