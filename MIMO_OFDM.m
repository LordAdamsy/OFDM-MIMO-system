% ͨ��ϵͳ������SOC���� 2023���＾ѧ��
clear;close all;clc;

% SNR��λΪdB�������е���
SNR_list = 0 : 3 : 30;

% ÿ��SNR���������OFDM���������������е���
block_num = 3000;

% FFT����
N_fft = 1024;     

% ���ƽ�����4(QPSK)��16(16QAM)
Q = 16;

% ÿ�����ų��صı�����
B = log2(Q);

% ����ؼ���
BER_count = zeros(size(SNR_list));

% ������Ϣ�����ز�����
N_sc = 896;

% ���䣬������������
N_t = 3;
N_r = 3;

% ÿ��OFDM����������ı�����
N_bit = N_sc * B * N_t;

% CP����
length_CP = 73;

% 0:ZF 1:MMSE
equal_method = 0;

% ��ȡ��Ч�����ŵ��弤��Ӧhij������i=1,2,3��j=1,2,3
% hij��ʾ�ӵ�j���������ߵ���i���������ߵĵ�Ч�����ŵ��弤��Ӧ
load('h.mat');

for snr_count = 1 : length(SNR_list)
    
    snr = SNR_list(snr_count);
    
    % ���￼�Ƿ�������֮��ĵȹ��ʷ��䣬��ÿ������������ÿ���������ز��Ϸ�����ŵĹ���Ϊ1 / Nt
    % ��������������и��������ߵ�Ƶ������ŵ�ƽ������֮��Ϊ1��SNR�Ķ���Ϊ��N_sc / (N_fft ^ 2 * ��������)
    sigma2 = N_sc / (10 ^ (snr / 10) * N_fft ^ 2);
    
    for count = 1 : block_num
        snr_count, count
        
        %������Ϣ����
        msg = round(rand(1, N_bit));
        
        %% MIMO-OFDM���Ʋ���
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % �������Ĵ��룬��ɴӱ��ص�3·����������������QPSK��16QAM��������ӳ��,��Ƶ����ɻ���SVD��Ԥ���룬����OFDM���Ʋ����CP
        % ���룺��Ϣ��������msg
        % �����3·OFDM���ƺ����CP��ʱ���ź�����x_ofdm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % qpsk/16qam modulation
        if(Q == 4)
            msg_mod = (1/sqrt(N_t))*pskmod(bit2int(reshape(msg, B, []), B), Q);
        elseif(Q == 16)
            msg_mod = (1/sqrt(N_t))*(1/sqrt(10))*qammod(bit2int(reshape(msg, B, []), B), Q);
        end
        
        % divide data into 3 data streams
        msg_3_tr = reshape(msg_mod, [], 3).';
        msg_3_tr_zeros = [msg_3_tr, zeros(N_t, N_fft-N_sc)];
        
        % channel matrix
        H = zeros(N_r, N_t, N_fft);
        H(1, 1, :) = fft(h11, N_fft);
        H(1, 2, :) = fft(h12, N_fft);
        H(1, 3, :) = fft(h13, N_fft);
        H(2, 1, :) = fft(h21, N_fft);
        H(2, 2, :) = fft(h22, N_fft);
        H(2, 3, :) = fft(h23, N_fft);
        H(3, 1, :) = fft(h31, N_fft);
        H(3, 2, :) = fft(h32, N_fft);
        H(3, 3, :) = fft(h33, N_fft);

        % precoding
        U = zeros(N_r, N_r, N_fft);
        S = zeros(N_r, N_r, N_fft);
        msg_prec = zeros(N_r, N_fft);

        for i = 1:N_fft
            [U(:, :, i), S(:, :, i), V] = svd(H(:, :, i));
            msg_prec(:, i) = V * msg_3_tr_zeros(:, i);
        end
        
        % OFDM + D to S
        msg_ofdm = ifft(msg_prec, [], 2);
        
        % cyclic prefix
        x_ofdm = [msg_ofdm(:, length(msg_ofdm)-length_CP+1:end), msg_ofdm]; 

        %% �ŵ����䲿��
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % �������Ĵ��룬�������x_ofdm�����ྶ�ŵ�hij������նˣ������AWGN��ȥ��CP��ע��ʵ�鲿�������ʸ�Ϊ�������ʵ�һ��
        % �����Ѿ������CP����������һ��OFDM���ŶԱ�OFDM���ŵ�Ӱ��
        % ���룺3·���CP���ʱ��OFDM�������x_ofdm
        % �����3·����ΪN_fft�Ľ���OFDM����r_ofdm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % r_1 = conv(x_ofdm(1, :), h11, 'same') + conv(x_ofdm(2, :), h12, 'same') + conv(x_ofdm(3, :), h13, 'same');
        % r_2 = conv(x_ofdm(1, :), h21, 'same') + conv(x_ofdm(2, :), h22, 'same') + conv(x_ofdm(3, :), h23, 'same');
        % r_3 = conv(x_ofdm(1, :), h31, 'same') + conv(x_ofdm(2, :), h32, 'same') + conv(x_ofdm(3, :), h33, 'same');
        
        r_1 = conv(x_ofdm(1, :), h11) + conv(x_ofdm(2, :), h12) + conv(x_ofdm(3, :), h13);
        r_2 = conv(x_ofdm(1, :), h21) + conv(x_ofdm(2, :), h22) + conv(x_ofdm(3, :), h23);
        r_3 = conv(x_ofdm(1, :), h31) + conv(x_ofdm(2, :), h32) + conv(x_ofdm(3, :), h33);

        r_time = [r_1; r_2; r_3];      

        r_channel = r_time + randn(N_r, length(r_1)) * sqrt(0.5 * sigma2) * (1+1i);
        r_ofdm = r_channel(:, length_CP+1:length_CP+N_fft);
        
        %% MIMO-OFDM�������       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ��������Լ��Ĵ��룬�Խ��յ���OFDM-MIMO���Ž��������MMSE������о������㲢�ָ��������
        % ���룺3·����OFDM����r_ofdm
        % ����������ı�������msg_r
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        r_fft = fft(r_ofdm, [], 2);

        r_com = zeros(N_r, N_fft);
        r_eq = zeros(N_r, N_fft);
        % combine
        for i = 1:N_fft
            r_com(:, i) = U(:, :, i)' * r_fft(:, i);        
            % equalization
            if(equal_method == 0)
                r_eq(:, i) = inv(S(:, :, i)) * r_com(:, i);
            elseif(equal_method == 1)
                r_eq(:, i) = inv((S(:, :, i)' * S(:, :, i) + sigma2 * eye(N_t))) * S(:, :, i)' * r_com(:, i);
            end
        end
        
        % delete zeros
        r_data = sqrt(N_t) * reshape(r_eq(:, 1:N_sc).', [], 1).';
        
        % demod
        if(Q == 4)
            msg_demod = pskdemod(r_data, Q);
        elseif(Q == 16)
            msg_demod = qamdemod(sqrt(10) * r_data, Q);
        end
        
        msg_r = reshape(int2bit(msg_demod, B), 1, []);
        
        %% �������ͳ��
        BER_count(snr_count) = BER_count(snr_count) + sum(abs(msg_r - msg));
        
    end
end

% �����ʼ���
BER = BER_count / (block_num * N_bit);

%% BER���߻���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������Ĵ��룬���ð�������꣬ʹ��semilogy��������BER-SNR����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(SNR_list, BER, 'LineWidth', 1.5);
xlabel('SNR/dB');
ylabel('BER'); 
grid on;

