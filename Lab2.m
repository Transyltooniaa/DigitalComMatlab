close all;
clear all;
clc;

blockLength = 1000;
nBlocks = 20000;
SNRdb = 0:1:15;              % Signal-to-noise ratio in dB
SNR = 10.^(SNRdb/10);        % Convert SNR from dB to linear scale
No = 1;                      % Noise power spectral density
BER_BPSK = zeros(1, length(SNRdb)); % Initialize BER array

for blk = 1:nBlocks
    Bits = randi([0,1],1,blockLength);  % 1. Generate random bits of length equal to block length
    Sym = 2*Bits - 1;                     % 2. Generate BPSK symbols (map bit 1 -> +1, bit 0 -> -1)
    noise = sqrt(No/2) * (randn(1,blockLength)+1j*randn(1,blockLength)); % 3. Generate noise with variance No/2
   
    for K = 1:length(SNRdb)
        TxSym = sqrt(SNR(K)) * Sym;         % 4. Transmit the symbols, adjust amplitude according to SNR
        RxSym = TxSym + noise;              % 5. Calculate the received symbols (add noise)
        DecBits = real((RxSym) > 0);                % 6. Decision rule: if RxSym > 0, decode as 1, else decode as 0
        BER_BPSK(K) = BER_BPSK(K) + sum(DecBits ~= Bits)/blockLength;  % 7. Calculate BER for each SNR value
    end
end

BER = BER_BPSK / nBlocks;  % Compute average BER

semilogy(SNRdb,BER,'g','linewidth',2.0,'MarkerSize',9.0);
grid on;
legend('BER');
xlabel('SNR(db)');
ylabel('BER')
