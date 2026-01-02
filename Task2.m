% Task 2: Fixed SNR = 10 dB, sweep Nr from 3 to 6, Nt = 3
Nt = 3;
SNR_fixed = 10;
 numBits =1e5;  
Nr_list = 3:6;

BER_pts_ML        = zeros(1,length(Nr_list));
BER_pts_ZF        = zeros(1,length(Nr_list));
BER_pts_ZF_SIC    = zeros(1,length(Nr_list));
BER_pts_MMSE      = zeros(1,length(Nr_list));
BER_pts_MMSE_SIC  = zeros(1,length(Nr_list));

for i = 1:length(Nr_list)
    Nr = Nr_list(i);
    BER_pts_ML(i)        = spatial_mux_ML(SNR_fixed, numBits, Nt, Nr);
    BER_pts_ZF(i)        = spatial_mux_zf(SNR_fixed, numBits, Nt, Nr);
    BER_pts_ZF_SIC(i)    = spatial_mux_sic(SNR_fixed, numBits, Nt, Nr);
    BER_pts_MMSE(i)      = spatial_mux_mmse(SNR_fixed, numBits, Nt, Nr);
    BER_pts_MMSE_SIC(i)  = spatial_mux_mmse_sic(SNR_fixed, numBits, Nt, Nr);
end

figure;
semilogy(Nr_list, BER_pts_ML,       'k-o','LineWidth',2); hold on;
semilogy(Nr_list, BER_pts_ZF,       'b-s','LineWidth',2);
semilogy(Nr_list, BER_pts_ZF_SIC,   'r-d','LineWidth',2);
semilogy(Nr_list, BER_pts_MMSE,     'm-^','LineWidth',2);
semilogy(Nr_list, BER_pts_MMSE_SIC, 'g-v','LineWidth',2);
grid on;
xlabel('Number of Receivers (Nr)');
ylabel('BER at SNR = 10 dB');
legend('ML','ZF','ZF-SIC','MMSE','MMSE-SIC','Location','southwest');
title('BER vs Nr (3Tx, SNR=10 dB) for ML, ZF, ZF-SIC, MMSE, MMSE-SIC');

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

