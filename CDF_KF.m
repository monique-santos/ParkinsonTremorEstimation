clc; clear; close all;

%% 1. CARREGAMENTO DOS DADOS
[file, path] = uigetfile('*.txt','Selecione os dados do Giroscopio');
if isequal(file,0); return; end
dados = load(fullfile(path,file));

t_orig = dados(:,1);
gyrox  = dados(:,5);
N  = length(t_orig);
Fs = 1/mean(diff(t_orig));
Ts = 1/Fs;
t = (0:Ts:(N-1)*Ts)';

%% =========================================================
% 2. MOVIMENTO VOLUNTÁRIO - CDF
% =========================================================
theta = 0.975;
g = 1 - theta^2;
h = (1 - theta)^2;
x_p = gyrox(1);
v_p = 0;
ref_vol = zeros(N,1);

for k = 1:N
    residuo = gyrox(k) - x_p;
    x_atual = x_p + g * residuo;
    v_atual = v_p + (h/Ts) * residuo;
    ref_vol(k) = x_atual;
    % Predição
    x_p = x_atual + Ts * v_atual;
    v_p = v_atual;
end

%% =========================================================
% 3. FILTRO BENEDICT-BORDNER
% =========================================================
alpha = 0.041;
g_bb = alpha;
h_bb = (alpha^2) / (2 - alpha);
x_bb = zeros(N,1);
v_p_bb = 0;
x_p_bb = gyrox(1);

for k = 1:N
    res_bb = gyrox(k) - x_p_bb;
    x_bb(k) = x_p_bb + g_bb * res_bb;
    v_atual_bb = v_p_bb + (h_bb/Ts) * res_bb;
    x_p_bb = x_bb(k) + Ts * v_atual_bb;
    v_p_bb = v_atual_bb;
end

%% =========================================================
% 4. FILTRO DE KALMAN (MOVIMENTO VOLUNTÁRIO)
% =========================================================
A_kf_v = [1 Ts; 0 1];
H_kf_v = [1 0];
Q_kf_v = [(Ts^4)/4 (Ts^3)/2; (Ts^3)/2 Ts^2] * 0.0001;
R_kf_v = 0.0001;

x_kf_vol = zeros(N,1); % RENOMEADO para evitar conflito
X_hat_kf_v = [gyrox(1); 0];
P_kf_v = eye(2);

for k = 1:N
    % Predição
    X_pred = A_kf_v * X_hat_kf_v;
    P_pred = A_kf_v * P_kf_v * A_kf_v' + Q_kf_v;
    % Ganho de Kalman
    K = (P_pred * H_kf_v') / (H_kf_v * P_pred * H_kf_v' + R_kf_v);
    % Atualização
    X_hat_kf_v = X_pred + K * (gyrox(k) - H_kf_v * X_pred);
    P_kf_v = (eye(2) - K * H_kf_v) * P_pred;
    x_kf_vol(k) = X_hat_kf_v(1);
end

%% =========================================================
% FILTRAGEM DE REFERÊNCIA
% =========================================================
fc = 2; 
[b,a] = butter(4, fc/(Fs/2));
gyrox_filtrado = filtfilt(b, a, gyrox);

%% =========================================================
% MÉTRICA KTE
% =========================================================
calc_kte = @(original, estimado) ...
    sqrt( mean(abs(original - estimado))^2 + ...
    var(abs(original - estimado))^2 );

kte_cdf = calc_kte(gyrox_filtrado, ref_vol);
kte_bb  = calc_kte(gyrox_filtrado, x_bb);
kte_kf  = calc_kte(gyrox_filtrado, x_kf_vol);

disp('========== KTE ==========')
fprintf('CDF: %.6f\n', kte_cdf);
fprintf('Benedict-Bordner: %.6f\n', kte_bb);
fprintf('Kalman: %.6f\n', kte_kf);

%% =========================================================
% 5. ESTIMAÇÃO DO TREMOR (WFLC + KALMAN EM CASCATA)
% =========================================================
tremor_bruto = gyrox - ref_vol;

% Configuração WFLC
mu0_wflc = 1e-4;      
mu1_wflc = 0.1;       
mub_wflc = 0.01;      
w0_est = 2 * pi * 6;  
W_wflc = zeros(2, 1);
bias_wflc = 0;

% Configuração KF para Tremor
x_kf_tremor = [w0_est; 0; 0; 0]; 
F_func = @(w) [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, cos(w), sin(w), 0];
H_kf_t = [0, 0, 0, 1];
R_t = 0.01;
Q_t = diag([1e-6, 1e-2, 1e-2, 1e-4]);
P_kf_t = eye(4) * 0.1;

tremor_estimado_kf = zeros(N, 1);

for k = 1:N
    % --- WFLC ---
    theta_k = w0_est * k * Ts;
    Xk_wflc = [sin(theta_k); cos(theta_k)];
    tremor_est_wflc = (W_wflc' * Xk_wflc) + bias_wflc;
    erro_wflc = tremor_bruto(k) - tremor_est_wflc;
    
    sum_term = (W_wflc(1) * Xk_wflc(2) - W_wflc(2) * Xk_wflc(1));
    w0_est = w0_est + 2 * mu0_wflc * erro_wflc * sum_term;
    w0_est = max(2*pi*3, min(2*pi*12, w0_est)); % Trava fisiológica
    
    W_wflc = W_wflc + 2 * mu1_wflc * erro_wflc * Xk_wflc;
    bias_wflc = bias_wflc + 2 * mub_wflc * erro_wflc;
    
    % --- Kalman Tremor ---
    F_k = F_func(w0_est * Ts);
    x_pred = F_k * x_kf_tremor;
    P_pred = F_k * P_kf_t * F_k' + Q_t;
    
    inovacao = tremor_bruto(k) - H_kf_t * x_pred;
    S_mat = H_kf_t * P_pred * H_kf_t' + R_t;
    K_gain = P_pred * H_kf_t' / S_mat;
    
    x_kf_tremor = x_pred + K_gain * inovacao;
    P_kf_t = (eye(4) - K_gain * H_kf_t) * P_pred;
    
    tremor_estimado_kf(k) = x_kf_tremor(4);
end

%% =========================================================
% 6. PLOTS NO TEMPO
% =========================================================
figure('Name','Movimento Bruto vs Estimativas','Color','w');
plot(t, gyrox, 'k', 'HandleVisibility','on'); hold on;
plot(t, ref_vol, 'b', 'LineWidth', 1.5);
plot(t, x_bb, 'g', 'LineWidth', 1.5);
plot(t, x_kf_vol, 'r', 'LineWidth', 1.5);
title('Estimativas do Movimento Voluntário');
ylabel('Velocidade Angular (rad/s)');
xlabel('Tempo (s)');
xlim([0 10]); ylim([-6 6]);
legend('Bruto', 'CDF', 'Benedict-Bordner', 'Kalman');
grid on;

figure('Name','Tremor Isolado','Color','w');
plot(t, tremor_bruto, 'Color',[0.7 0.7 0.7]); hold on;
plot(t, tremor_estimado_kf, 'r', 'LineWidth',1.5);
title('Tremor Estimado (WFLC + KF)');
xlabel('Tempo (s)');
ylabel('Velocidade Angular (rad/s)');
legend('Resíduo Bruto', 'Tremor Filtrado');
grid on;

%% =========================================================
% 7. ANÁLISE DE FREQUÊNCIA (FFT)
% =========================================================
calc_fft = @(s) abs(fft(s)/N);
f_axis = (0:N-1)*(Fs/N);
cutoff = round(N/2);

fft_bruto = calc_fft(gyrox);
fft_vol_cdf = calc_fft(ref_vol);
fft_vol_kf = calc_fft(x_kf_vol);

figure('Name','Análise de Frequência','Color','w');
plot(f_axis(1:cutoff), fft_bruto(1:cutoff), 'k--', 'LineWidth',1); hold on;
plot(f_axis(1:cutoff), fft_vol_cdf(1:cutoff), 'b', 'LineWidth',1.5);
plot(f_axis(1:cutoff), fft_vol_kf(1:cutoff), 'r', 'LineWidth',1.5);
title('Espectro de Magnitude');
xlabel('Frequência (Hz)');
ylabel('|X(f)|');
legend('Sinal Bruto','Estimativa CDF', 'Estimativa Kalman');
xlim([0 15]);
ylim([0 0.6]);
grid on;

figure('Name','Sinal Bruto Individual','Color','w');
plot(t, gyrox, 'k');
title('Sinal Bruto Adquirido');
xlabel('Tempo (s)'); ylabel('rad/s');
xlim([0 10]); ylim([-6 6]);
grid on;

figure('Name','FFT do Sinal Bruto','Color','w');
plot(f_axis(1:cutoff), fft_bruto(1:cutoff), 'k', 'LineWidth', 1.5);
title('Espectro de Magnitude - Sinal Bruto');
xlabel('Frequência (Hz)');
ylabel('|X(f)|');
xlim([0 15]);
grid on;

