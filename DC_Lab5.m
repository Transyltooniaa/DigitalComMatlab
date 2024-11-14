close all;
clear all;
clc;

% Parameters
blocklength = 1000;  % Block length
nBlocks = 5000;     % Number of blocks
SNRdb = 0:1:15;      % SNR in dB
SNR = 10.^(SNRdb/10); % Linear scale SNR
No = 1;              % Noise power spectral density (No = 1 for simplicity)
L = 5;               % Number of receive antennas
BER_BPSK = zeros(1, length(SNRdb));  % Initialize BER array

% Calculate Theoretical BER
for K = 1:length(SNRdb)
    BER_theoretical(K) = 0.5 * qfunc(sqrt(SNR(K))); % Calculate theoretical BER
end

% Simulation
for blk = 1:nBlocks
    % 1. Generate random bits of length equal to block length
    Bits = randi([0 1], 1, blocklength);  
    
    % 2. Generate BPSK symbols: map bit 1----> +1 and bit 0----> -1
    sym = 2 * Bits - 1; 
    
    % Generate noise for each SNR value
    noise = sqrt(No/2) * (randn(L, blocklength) + 1j * randn(L, blocklength)); % L x blocklength noise matrix
    
    % Generate Rayleigh fading channel for L receive antennas
    h = (randn(L, blocklength) + 1j * randn(L, blocklength)) / sqrt(2); % Rayleigh channel matrix
    
    for K = 1:length(SNRdb)
        % 3. Transmit the symbols, adjust amplitude relative to desired SNR
        TxSym = sqrt(SNR(K)) * sym;
        
        % 4. Calculate the received symbols (h * TxSym + noise)
        RxSym = h .* TxSym + noise; % L x blocklength received signals
        
        % 5. MRC Combining: sum the received signals from all antennas
        combined_signal = sum(conj(h) .* RxSym, 1); % Sum across antennas (L dimension)
        
        % 6. Make decision based on the combined signal to recover the bits
        DecBits = real(combined_signal) > 0;  % If real part > 0, bit is '1', else bit is '0'
        
        % 7. Calculate BER for each iteration
        error_count = sum(Bits ~= DecBits); % Total errors for the current block
        
        % Update BER
        BER_BPSK(K) = BER_BPSK(K) + error_count;
    end
end

% Compute Average BER
BER_BPSK = BER_BPSK / (blocklength * nBlocks);  % Average over blocks and length

% Plot the BER vs SNR graph
figure;
semilogy(SNRdb, BER_BPSK, 'g-o', 'linewidth', 2.0, 'MarkerSize', 9.0, 'DisplayName', 'Simulated BER');
hold on;
semilogy(SNRdb, BER_theoretical, 'b--', 'linewidth', 2.0, 'DisplayName', 'Theoretical BER');
grid on;
legend('show');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for BPSK over Rayleigh Fading Channel with MRC (L = 5)');
grid on;

function y = qfunc(x)
    % Q-function definition
    y = 0.5 * erfc(x / sqrt(2));
end
