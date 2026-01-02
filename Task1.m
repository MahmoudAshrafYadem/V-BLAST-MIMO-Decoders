
% Task 1: 3x3 MIMO, BER vs SNR
Nt = 3;
Nr = 3;
numBits = 1e5;
SNRdB = 0:2:20;

BER_ML    = zeros(1,length(SNRdB));
BER_ZF    = zeros(1,length(SNRdB));
BER_ZF_SIC= zeros(1,length(SNRdB));
BER_MMSE  = zeros(1,length(SNRdB));
BER_MMSE_SIC = zeros(1,length(SNRdB));

for k = 1:length(SNRdB)
    BER_ML(k)        = spatial_mux_ML(SNRdB(k), numBits, Nt, Nr);
    BER_ZF(k)        = spatial_mux_zf(SNRdB(k), numBits, Nt, Nr);
    BER_ZF_SIC(k)    = spatial_mux_sic(SNRdB(k), numBits, Nt, Nr);
    BER_MMSE(k)      = spatial_mux_mmse(SNRdB(k), numBits, Nt, Nr);
    BER_MMSE_SIC(k)  = spatial_mux_mmse_sic(SNRdB(k), numBits, Nt, Nr);
end

figure;
semilogy(SNRdB, BER_ML,       'k-o','LineWidth',2); hold on;
semilogy(SNRdB, BER_ZF,       'b-s','LineWidth',2);
semilogy(SNRdB, BER_ZF_SIC,   'r-d','LineWidth',2);
semilogy(SNRdB, BER_MMSE,     'm-^','LineWidth',2);
semilogy(SNRdB, BER_MMSE_SIC, 'g-v','LineWidth',2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('ML','ZF','ZF-SIC','MMSE','MMSE-SIC','Location','southwest');
title('BER vs SNR for 3x3 MIMO Spatial Multiplexing (BPSK)');

function BER_spatial_ML = spatial_mux_ML(SNR_dB, numBits, Nt, Nr)

    % Generate all possible BPSK symbol vectors for Nt streams
    bpsk = [-1 1];
    M = 2^Nt; % number of candidates
    X = zeros(Nt, M);
    for m = 0:M-1
        bits_m = dec2bin(m,Nt) - '0';   % binary vector length Nt
        X(:,m+1) = 2*bits_m.' - 1;      % map {0,1} -> {-1,+1}
    end

    SNR = 10^(SNR_dB/10);
    noise_var = 1/(2*SNR);
    bit_errors = 0;

    for k = 1:numBits
        % Random transmit bits
        bits = randi([0 1], Nt, 1);
        x = 2*bits - 1;

        % Channel and noise
        H = (randn(Nr,Nt) + 1j*randn(Nr,Nt))/sqrt(2);
        n = sqrt(noise_var)*(randn(Nr,1) + 1j*randn(Nr,1));
        y = H*x + n;

        % Compute received candidates
        Y_hat = H * X;   % size Nr × M
        Y = repmat(y, 1, M);

        % ML decision: minimum distance
        distances = sum(abs(Y - Y_hat).^2, 1);
        [~, idx] = min(distances);
        x_hat = X(:, idx);

        bits_hat = (x_hat + 1)/2;
        bit_errors = bit_errors + sum(bits ~= bits_hat);
    end

    BER_spatial_ML = bit_errors / (numBits * Nt);
end


function BER_spatial_zf = spatial_mux_zf(SNR_dB, numBits, Nt, Nr)

    SNR = 10^(SNR_dB/10);
    noise_var = 1/(2*SNR);
    bit_errors = 0;

    for k = 1:numBits

        bits = randi([0 1], Nt, 1);     
        x = 2*bits - 1;               

        H = (randn(Nr,Nt) + 1j*randn(Nr,Nt))/sqrt(2);
        n = sqrt(noise_var)*(randn(Nr,1) + 1j*randn(Nr,1));

        y = H*x + n;

        W_zf = (H'*H)\H';       

        x_hat = W_zf * y;

        bits_hat = real(x_hat) > 0;
        bit_errors = bit_errors + sum(bits ~= bits_hat);

    end

    BER_spatial_zf = bit_errors / (numBits * Nt);
end



function BER_spatial_sic = spatial_mux_sic(SNR_dB, numBits, Nt, Nr)

    SNR = 10^(SNR_dB/10);
    noise_var = 1/(2*SNR);
    bit_errors = 0;

    for k = 1:numBits

    
        bits_original = randi([0 1], Nt, 1);
        x = 2*bits_original - 1;          

        
        H = (randn(Nr,Nt) + 1j*randn(Nr,Nt))/sqrt(2);
        n = sqrt(noise_var)*(randn(Nr,1) + 1j*randn(Nr,1));

       
        y = H*x + n;

        y_sic = y;                
        H_sic = H;                
        bits_hat = zeros(Nt,1);   
        stream_indices = 1:Nt;    

        for m = 1:Nt

            W_zf = (H_sic'*H_sic)\H_sic';

            
            noise_enhancement = diag(W_zf*W_zf');
            [~, idx_min] = min(noise_enhancement);  

            x_hat = W_zf(idx_min,:) * y_sic;
            bits_hat(stream_indices(idx_min)) = real(x_hat) > 0;

            y_sic = y_sic - H_sic(:,idx_min)*(2*bits_hat(stream_indices(idx_min))-1);

            H_sic(:,idx_min) = [];
            stream_indices(idx_min) = [];
        end

        bit_errors = bit_errors + sum(bits_hat ~= bits_original);

    end

    BER_spatial_sic = bit_errors / (numBits * Nt);
end
function BER_spatial_mmse = spatial_mux_mmse(SNR_dB, numBits, Nt, Nr)

    SNR = 10^(SNR_dB/10);
    noise_var = 1/(2*SNR);
    bit_errors = 0;

    for k = 1:numBits

        bits = randi([0 1], Nt, 1);
        x = 2*bits - 1;

        H = (randn(Nr,Nt) + 1j*randn(Nr,Nt))/sqrt(2);
        n = sqrt(noise_var)*(randn(Nr,1) + 1j*randn(Nr,1));
        y = H*x + n;

        % MMSE linear filter
        W_mmse = (H'*H + noise_var*eye(Nt)) \ H';
        x_hat = W_mmse * y;

        bits_hat = real(x_hat) > 0;
        bit_errors = bit_errors + sum(bits ~= bits_hat);

    end

    BER_spatial_mmse = bit_errors / (numBits * Nt);
end


function BER_spatial_mmse_sic = spatial_mux_mmse_sic(SNR_dB, numBits, Nt, Nr)

    SNR = 10^(SNR_dB/10);
    noise_var = 1/(2*SNR);
    bit_errors = 0;

    for k = 1:numBits

        bits_original = randi([0 1], Nt, 1);
        x = 2*bits_original - 1;

        H = (randn(Nr,Nt) + 1j*randn(Nr,Nt))/sqrt(2);
        n = sqrt(noise_var)*(randn(Nr,1) + 1j*randn(Nr,1));
        y = H*x + n;

        y_sic = y;
        H_sic = H;
        bits_hat = zeros(Nt,1);
        stream_indices = 1:Nt;

        for m = 1:Nt
            % MMSE filter for current sub-channel
            W_mmse = (H_sic'*H_sic + noise_var*eye(size(H_sic,2))) \ H_sic';

            % Choose stream with minimum post-filter noise enhancement (or max SINR proxy)
            enh = diag(W_mmse*W_mmse');      % proxy for noise enhancement
            [~, idx_best] = min(enh);

            % Detect chosen stream
            x_hat_m = W_mmse(idx_best,:) * y_sic;
            detected_bit = real(x_hat_m) > 0;
            bits_hat(stream_indices(idx_best)) = detected_bit;

            % Cancel its contribution from y
            x_cancel = 2*detected_bit - 1;
            y_sic = y_sic - H_sic(:,idx_best) * x_cancel;

            % Remove the decoded column
            H_sic(:,idx_best) = [];
            stream_indices(idx_best) = [];
        end

        bit_errors = bit_errors + sum(bits_hat ~= bits_original);

    end

    BER_spatial_mmse_sic = bit_errors / (numBits * Nt);
end

