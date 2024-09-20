close all;
clear all;
clc;

% Parameters
blocklength = 1000;  % Block length
nBlocks = 20000;     % Number of blocks
SNRdb = 0:1:15;      % SNR in dB
SNR = 10.^(SNRdb/10); % Linear scale SNR
No = 1;              % Noise power spectral density (No = 1 for simplicity)
BER_BPSK = zeros(1, length(SNRdb));  % Initialize BER array
BER_theoretical = zeros(1, length(SNRdb));  % Initialize theoretical BER array

% Calculate Theoretical BER
for K = 1:length(SNRdb)
    BER_theoretical(K) = 0.5 * qfunc(sqrt(SNR(K))); % Calculate theoretical BER
end

% Simulation
for blk = 1:nBlocks
    % 1. Generate random bits of length equal to block length
    Bits = randi([0 1], 1, blocklength);  
    
    % 2. Generate BPSK symbols: map bit 1----> +1 and bit 0----> -1
    sym = 2*Bits - 1; 
    
    % Generate noise for each SNR value
    noise = sqrt(No/2) * (randn(1, blocklength) + 1j*randn(1, blocklength));
    
    % Generate Rayleigh fading channel
    h = (randn(1, blocklength) + 1j*randn(1, blocklength)) / sqrt(2); % Rayleigh channel
    
    for K = 1:length(SNRdb)
        % 3. Transmit the symbols, adjust amplitude relative to desired SNR
        TxSym = sqrt(SNR(K)) * sym;
        
        % 4. Calculate the received symbols (h * TxSym + noise)
        RxSym = h .* TxSym + noise;
        
        % 5. Equalize the received symbols (i.e., divide by the channel gain h)
        EqualizedSym = RxSym ./ h;
        
        % 6. Make decision based on the equalized symbols to recover the bits
        DecBits = real(EqualizedSym) > 0;  % If real part > 0, bit is '1', else bit is '0'
        
        % 7. Calculate BER for each iteration
        BER_BPSK(K) = BER_BPSK(K) + sum(Bits ~= DecBits);
    end
end

% Compute Average BER
BER_BPSK = BER_BPSK / (blocklength * nBlocks);  % Average over blocks and length

% Plot the BER vs SNR graph
figure;
semilogy(SNRdb, BER_BPSK, 'g-o', 'linewidth', 2.0, 'MarkerSize', 9.0);
hold on;
semilogy(SNRdb, BER_theoretical, 'b--', 'linewidth', 2.0);
grid on;
legend('Simulated BER', 'Theoretical BER');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for BPSK over Rayleigh Fading Channel');
grid on;

function y = qfunc(x)
    % Q-function definition
    y = 0.5 * erfc(x / sqrt(2));
end
