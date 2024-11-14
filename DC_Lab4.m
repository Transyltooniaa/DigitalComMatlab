clc;
clear;

nPilot = 10;
nData = 1e3;
blockLength = nPilot + nData;
nBlocks = 2e3;
SNRdb = 0:0.5:15;
SNR = 10.^(SNRdb / 10);


No = 1;  % Noise power
BER_pCSI = zeros(size(SNRdb)); % BER of perfect channel
BER_imCSI = zeros(size(SNRdb)); % BER of imperfect or estimated channel
hEstErr = zeros(size(SNRdb)); % Channel estimation error

for blk = 1:nBlocks
    h = (randn + 1i*randn) / sqrt(2);  % 1. Generate complex Gaussian channel (Rayleigh fading)
    Bits = randi([0 1], blockLength, 1); % 2. Generate bits of length equal to block length
    Sym = 2*Bits - 1;  % 3. BPSK modulation (mapping 0 -> -1, 1 -> +1)
    pSym = Sym(1:nPilot);  % 4. Get nPilot number of pilot symbols
    pNorm = norm(pSym)^2;  % 5. Calculate the norm square of pilot symbols
    
    noise = (randn(blockLength, 1) + 1i*randn(blockLength, 1)) / sqrt(2);  % Generate complex Gaussian noise with variance N0
    
    for K = 1:length(SNRdb)
        % Transmit symbols with noise
        TxSym = sqrt(SNR(K)) * Sym;  % Scale symbols by sqrt(SNR)
        RxSym = h * TxSym + sqrt(No) * noise;  % Received symbols (h is the channel, plus noise)
        
        % 6. Get the received observation vector corresponding to the pilot symbols
        pObser = RxSym(1:nPilot);
        
        % 7. Estimate channel using least squares (LS) estimator
        hEst = (pSym' * pObser) / pNorm;  
        
        % 8. Calculate channel estimation error
        hEstErr(K) = hEstErr(K) + abs(h - hEst)^2;
        
        % 9. Equalization for perfect CSI (perfect channel knowledge)
        EqSym_pCSI = RxSym / h;
        Decbits_pCSI = real(EqSym_pCSI) > 0;  % Decision rule (BPSK demodulation)
        BER_pCSI(K) = BER_pCSI(K) + sum(Decbits_pCSI(nPilot+1:end) ~= Bits(nPilot+1:end));
        
        % 10. Equalization for imperfect CSI (estimated channel knowledge)
        EqSym_imCSI = RxSym / hEst;
        Decbits_imCSI = real(EqSym_imCSI) > 0;  % Decision rule (BPSK demodulation)
        BER_imCSI(K) = BER_imCSI(K) + sum(Decbits_imCSI(nPilot+1:end) ~= Bits(nPilot+1:end));
    end
end

hEstMSE = hEstErr / nBlocks;  % Average channel estimation error (MSE)
BER_pCSI = BER_pCSI / (nData * nBlocks);  % Average BER for perfect CSI
BER_imCSI = BER_imCSI / (nData * nBlocks);  % Average BER for imperfect CSI

% Plotting results
figure;
plot(SNRdb, 10*log10(hEstMSE), 'b', 'LineWidth', 2.0)
grid on
legend('Channel Estimation MSE')
xlabel('SNR (dB)')
ylabel('MSE (dB)')

figure;
semilogy(SNRdb, BER_pCSI, 'b', 'LineWidth', 2.0)
hold on
semilogy(SNRdb, BER_imCSI, 'r', 'LineWidth', 2.0)
grid on
legend('BER with perfect CSI', 'BER with Estimated CSI')
xlabel('SNR (dB)')
ylabel('BER')
