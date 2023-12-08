% 通信系统仿真与SOC集成 2023年秋季学期
clear;close all;clc;

% SNR单位为dB，可自行调整
SNR_list = 0 : 3 : 30;

% 每个SNR下所传输的OFDM符号数量，可自行调整
block_num = 3000;

% FFT点数
N_fft = 1024;     

% 调制阶数：4(QPSK)，16(16QAM)
Q = 16;

% 每个符号承载的比特数
B = log2(Q);

% 误比特计数
BER_count = zeros(size(SNR_list));

% 传输信息的子载波数量
N_sc = 896;

% 发射，接收天线数量
N_t = 3;
N_r = 3;

% 每个OFDM符号所传输的比特数
N_bit = N_sc * B * N_t;

% CP点数
length_CP = 73;

% 0:ZF 1:MMSE
equal_method = 0;

% 读取等效基带信道冲激响应hij，其中i=1,2,3，j=1,2,3
% hij表示从第j个发射天线到第i个接收天线的等效基带信道冲激响应
load('h.mat');

for snr_count = 1 : length(SNR_list)
    
    snr = SNR_list(snr_count);
    
    % 这里考虑发射天线之间的等功率分配，即每个发射天线在每个非零子载波上发射符号的功率为1 / Nt
    % 计算噪声方差，其中各发射天线的频域发射符号的平均功率之和为1，SNR的定义为：N_sc / (N_fft ^ 2 * 噪声方差)
    sigma2 = N_sc / (10 ^ (snr / 10) * N_fft ^ 2);
    
    for count = 1 : block_num
        snr_count, count
        
        %生成信息比特
        msg = round(rand(1, N_bit));
        
        %% MIMO-OFDM调制部分
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 请添加你的代码，完成从比特到3路发射天线数据流的QPSK（16QAM）星座点映射,在频域完成基于SVD的预编码，进行OFDM调制并添加CP
        % 输入：信息比特序列msg
        % 输出：3路OFDM调制后并添加CP的时域信号序列x_ofdm
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

        %% 信道传输部分
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 请添加你的代码，发射符号x_ofdm经过多径信道hij到达接收端，并添加AWGN、去掉CP。注意实虚部噪声功率各为噪声功率的一半
        % 由于已经添加了CP，不考虑上一个OFDM符号对本OFDM符号的影响
        % 输入：3路添加CP后的时域OFDM发射符号x_ofdm
        % 输出：3路长度为N_fft的接收OFDM符号r_ofdm
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
        
        %% MIMO-OFDM解调部分       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 请添加你自己的代码，对接收到的OFDM-MIMO符号进行迫零或MMSE解调，判决星座点并恢复传输比特
        % 输入：3路接收OFDM符号r_ofdm
        % 输出：解调后的比特序列msg_r
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
        
        %% 误比特数统计
        BER_count(snr_count) = BER_count(snr_count) + sum(abs(msg_r - msg));
        
    end
end

% 误码率计算
BER = BER_count / (block_num * N_bit);

%% BER曲线绘制
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 请添加你的代码，采用半对数坐标，使用semilogy函数绘制BER-SNR曲线
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(SNR_list, BER, 'LineWidth', 1.5);
xlabel('SNR/dB');
ylabel('BER'); 
grid on;

