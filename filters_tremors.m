clc; clear; close all;

%% CARREGAMENTO E CONFIGURAÇÕES INICIAIS
[file, path] = uigetfile('*.txt');
if isequal(file,0)
    return;
end
dados = load(fullfile(path, file));

t_orig = dados(:,1);    
gyrox  = dados(:,5);    % Sinal do giroscópio (eixo X)
N  = length(t_orig);        
Ts = mean(diff(t_orig));    
Fs = 1 / Ts;           
t  = (0:Ts:(N-1)*Ts)';

fprintf('Frequência de amostragem: %.2f Hz\n', Fs);
fprintf('Número de amostras: %d\n', N);
fprintf('Duração do sinal: %.2f segundos\n', t(end));

%% FILTRO CDF - Remoção do Movimento Voluntário
theta = 0.975;
g = 1 - theta^2;
h = (1 - theta)^2;

x_cdf = zeros(N,1);     
v_cdf = zeros(N,1);     
x_pred = gyrox(1);
v_pred = 0;

for k = 1:N
    z = gyrox(k);
    x_cdf(k) = x_pred + g * (z - x_pred);
    v_cdf(k) = v_pred + (h / Ts) * (z - x_pred);
    x_pred = x_cdf(k) + Ts * v_cdf(k);
    v_pred = v_cdf(k);
end

tremor_entrada = gyrox - x_cdf;

%% WFLC - Weighted Frequency Fourier Linear Combiner 

M = 1;
mu0 = 1e-4;
mu1 = 0.1;
mub = 0.01;

f_min = 3.5;
w_min = 2 * pi * f_min;
w0_wflc = 2 * pi * 6;

W_wflc = zeros(2*M, 1);
bias_wflc = 0;
tremor_wflc = zeros(N,1);
freq_wflc = zeros(N,1);

for k = 1:N
    theta_k = w0_wflc * k * Ts;
    Xk = [sin(theta_k); cos(theta_k)];
    
    est_k = (W_wflc' * Xk) + bias_wflc;
    ek = tremor_entrada(k) - est_k;
    
    sum_term = (W_wflc(1) * Xk(2) - W_wflc(2) * Xk(1));
    w0_wflc = w0_wflc + 2 * mu0 * ek * sum_term;
    w0_wflc = max(w0_wflc, w_min);
    
    W_wflc = W_wflc + 2 * mu1 * ek * Xk;
    bias_wflc = bias_wflc + 2 * mub * ek;
    
    tremor_wflc(k) = est_k;
    freq_wflc(k) = w0_wflc / (2*pi);
end

%% BMFLC - Band-limited Multiple Fourier Linear Combiner 

M_bmflc = 1;
mu_bmflc = 0.01;
mub_bmflc = 0.001;
omega_0 = 2*pi*3;
omega_f = 2*pi*12;
G = 20;

delta_omega = (omega_f - omega_0) / (G + 1);
frequencias = omega_0 + (0:G-1) * delta_omega;

n_coef = 2 * M_bmflc * G;
W_bmflc = zeros(n_coef, 1);
bias_bmflc = 0;
tremor_bmflc = zeros(N, 1);
amplitudes = zeros(G, N);

for k = 1:N
    Xk = zeros(n_coef, 1);
    idx = 1;
    for r = 1:G
        omega_r = frequencias(r);
        theta_k = omega_r * k * Ts;
        
        for m = 1:M_bmflc
            Xk(idx) = sin(m * theta_k);
            idx = idx + 1;
            Xk(idx) = cos(m * theta_k);
            idx = idx + 1;
        end
    end
    
    estimativa = (W_bmflc' * Xk) + bias_bmflc;
    erro = tremor_entrada(k) - estimativa;
    
    W_bmflc = W_bmflc + 2 * mu_bmflc * erro * Xk;
    bias_bmflc = bias_bmflc + 2 * mub_bmflc * erro;
    
    tremor_bmflc(k) = estimativa;
    
    idx2 = 1;
    for r = 1:G
        a_r = W_bmflc(idx2);
        b_r = W_bmflc(idx2+1);
        amplitudes(r, k) = a_r^2 + b_r^2;
        idx2 = idx2 + 2;
    end
end

freq_bmflc = zeros(N, 1);
for k = 1:N
    numerador = 0;
    denominador = 0;
    for r = 1:G
        peso_amp = amplitudes(r, k);
        numerador = numerador + peso_amp * (frequencias(r)/(2*pi));
        denominador = denominador + peso_amp;
    end
    if denominador > 0
        freq_bmflc(k) = numerador / denominador;
    else
        freq_bmflc(k) = mean(frequencias)/(2*pi);
    end
end

%% KF - Filtro de Kalman em Cascata com WFLC 

M_kf = 1;
mu0_kf = 1e-4;
mu1_kf = 0.1;
mub_kf = 0.01;

w0_est_kf = 2 * pi * 6;
W_kf = zeros(2*M_kf, 1);
bias_kf = 0;

x_est = [w0_est_kf; 0; 0; 0];

F_kf = @(w_Ts) [1, 0, 0, 0;
                 0, 1, 0, 0;
                 0, 0, 1, 0;
                 0, cos(w_Ts), sin(w_Ts), 0];

H_kf = [0, 0, 0, 1];

R_kf = 0.01;
Q_kf = diag([1e-6, 1e-2, 1e-2, 1e-4]);

P_kf = eye(4) * 0.1;

tremor_kf = zeros(N, 1);
freq_kf = zeros(N, 1);

for k = 1:N
    theta_k = w0_est_kf * k * Ts;
    Xk_wflc = [sin(theta_k); cos(theta_k)];
    
    tremor_est_wflc = (W_kf' * Xk_wflc) + bias_kf;
    erro_wflc = tremor_entrada(k) - tremor_est_wflc;
    
    sum_term = (W_kf(1) * Xk_wflc(2) - W_kf(2) * Xk_wflc(1));
    w0_est_kf = w0_est_kf + 2 * mu0_kf * erro_wflc * sum_term;
    w0_est_kf = max(2*pi*3, min(2*pi*12, w0_est_kf));
    
    W_kf = W_kf + 2 * mu1_kf * erro_wflc * Xk_wflc;
    bias_kf = bias_kf + 2 * mub_kf * erro_wflc;
    
    freq_kf(k) = w0_est_kf / (2*pi);
    
    F_current = F_kf(w0_est_kf * Ts);
    x_pred = F_current * x_est;
    P_pred = F_current * P_kf * F_current' + Q_kf;
    
    inovacao = tremor_entrada(k) - H_kf * x_pred;
    S = H_kf * P_pred * H_kf' + R_kf;
    K_gain = P_pred * H_kf' / S;
    
    x_est = x_pred + K_gain * inovacao;
    P_kf = (eye(4) - K_gain * H_kf) * P_pred;
    
    tremor_kf(k) = x_est(4);
end

%% CÁLCULO DO FMSEd (Filtered Mean Square Error with delay correction)


function [fmse, delay_estimates, aligned_signal] = calcular_FMSEd(referencia, estimado, max_delay, mu_delay)
    % Parâmetros:
    %   referencia: sinal de referência (ground truth - tremor real)
    %   estimado: sinal estimado pelo algoritmo
    %   max_delay: atraso máximo a ser considerado (em amostras)
    %   mu_delay: passo de adaptação para estimação do atraso (LMS)
    
    if nargin < 3
        max_delay = round(length(referencia) * 0.1); 
    end
    if nargin < 4
        mu_delay = 0.01; % Passo de adaptação
    end
    
    N = length(referencia);
    
    % Inicialização do estimador de atraso adaptativo
    delay_estimates = zeros(N, 1);
    delay_estimates(1) = 0;
    
    % Filtro de atraso variável 
    aligned_signal = zeros(N, 1);
    
    for k = 1:N
        % Estimativa do atraso atual (baseado no erro quadrático)
        if k > 1
            % Algoritmo LMS para estimação do atraso
            % Derivada do erro em relação ao atraso
            if k > 2 && delay_estimates(k-1) > 0
                % Estimativa da diferença do sinal para cálculo do gradiente
                delta = (estimado(k-1) - estimado(min(k, N-1))) / Ts;
                % Atualização do atraso (minimização do erro quadrático)
                erro_instantaneo = referencia(k) - aligned_signal(k-1);
                delay_update = mu_delay * erro_instantaneo * delta;
                delay_estimates(k) = delay_estimates(k-1) + delay_update;
            else
                delay_estimates(k) = delay_estimates(k-1);
            end
            
            % Limitação do atraso
            delay_estimates(k) = max(0, min(max_delay, delay_estimates(k)));
        end
        
        % Aplicação do atraso compensado ao sinal estimado
        if delay_estimates(k) >= 1
            % Interpolação linear para atraso fracionário
            d_int = floor(delay_estimates(k));
            d_frac = delay_estimates(k) - d_int;
            
            if k - d_int >= 1
                if k - d_int - 1 >= 1
                    aligned_signal(k) = (1 - d_frac) * estimado(k - d_int) + ...
                                        d_frac * estimado(k - d_int - 1);
                else
                    aligned_signal(k) = estimado(k - d_int);
                end
            else
                aligned_signal(k) = estimado(k);
            end
        else
            % Interpolação para atraso sub-amostral
            if k > 1
                d_frac = delay_estimates(k);
                aligned_signal(k) = (1 - d_frac) * estimado(k) + d_frac * estimado(k-1);
            else
                aligned_signal(k) = estimado(k);
            end
        end
    end
    
    % Cálculo do FMSE 
    % FMSEd_k = sqrt(E[(s_k - t_{k - d_k})^2])
    erro_compensado = referencia - aligned_signal;
    fmse = sqrt(mean(erro_compensado.^2));
    
end

%% CÁLCULO DAS MÉTRICAS DE DESEMPENHO (Incluindo FMSEd)

fprintf('\n========== MÉTRICAS DE DESEMPENHO ==========\n');
  
fprintf('\n--- FMSEd (Filtered MSE com compensação de atraso) ---\n');



max_delay_samples = round(0.05 * Fs); 
mu_delay = 0.001; 

% Cálculo do FMSEd para cada algoritmo
[fmse_wflc, delay_wflc, aligned_wflc] = calcular_FMSEd(tremor_entrada, tremor_wflc, max_delay_samples, mu_delay);
[fmse_bmflc, delay_bmflc, aligned_bmflc] = calcular_FMSEd(tremor_entrada, tremor_bmflc, max_delay_samples, mu_delay);
[fmse_kf, delay_kf, aligned_kf] = calcular_FMSEd(tremor_entrada, tremor_kf, max_delay_samples, mu_delay);

fprintf('\n--- Resultados FMSEd ---\n');
fprintf('FMSEd WFLC:    %.6f rad/s (atraso médio: %.2f ms)\n', fmse_wflc, mean(delay_wflc)*1000/Fs);
fprintf('FMSEd BMFLC:   %.6f rad/s (atraso médio: %.2f ms)\n', fmse_bmflc, mean(delay_bmflc)*1000/Fs);
fprintf('FMSEd KF:      %.6f rad/s (atraso médio: %.2f ms)\n', fmse_kf, mean(delay_kf)*1000/Fs);

% Interpretação dos resultados
fprintf('\n--- Interpretação das Métricas ---\n');
fprintf(' FMSEd:\n');
fprintf('  WFLC: FMSEd=%.6f ', ...
     fmse_wflc);
fprintf('  BMFLC: FMSEd=%.6f ', ...
    fmse_bmflc);
fprintf('  KF:     FMSEd=%.6f ', ...
     fmse_kf);

% Frequência Dominante do Tremor
fft_residual = abs(fft(tremor_entrada)) / N;
f_axis = (0:floor(N/2)-1) * (Fs/N);
half_idx = 1:floor(N/2);
mag_residual = 2 * fft_residual(half_idx);

[~, idx_peak] = max(mag_residual(f_axis >= 3 & f_axis <= 12));
freq_dom = f_axis(find(f_axis >= 3 & f_axis <= 12, 1) + idx_peak - 1);
fprintf('\nFrequência Dominante do Tremor: %.2f Hz\n', freq_dom);

fprintf('\n=============================================\n');

%% PLOTS 

% Sinal Original vs Movimento Voluntário vs Tremor
figure('Name', 'Sinais no Domínio do Tempo', 'Color', 'w', 'Position', [50 50 1200 800]);

subplot(2,1,1);
plot(t, gyrox, 'k', 'LineWidth', 1);
title('Sinal Bruto do Giroscópio (Eixo X)');
ylabel('Velocidade Angular (rad/s)');
grid on; grid on; axis tight; 

subplot(2,1,2);
plot(t, tremor_entrada, 'k--', 'LineWidth', 0.8, 'DisplayName', 'Tremor bruto');
hold on;
plot(t, tremor_wflc, 'r', 'LineWidth', 1.2, 'DisplayName', 'WFLC');
plot(t, tremor_bmflc, 'b', 'LineWidth', 1.2, 'DisplayName', 'BMFLC');
plot(t, tremor_kf, 'g', 'LineWidth', 1.2, 'DisplayName', 'KF em Cascata');
title('Estimação do Tremor Patológico');
xlabel('Tempo (s)');
ylabel('Velocidade Angular (rad/s)');
legend('Location', 'best');
grid on; axis tight;
sgtitle('Análise Temporal dos Sinais');

% FFTs Comparativas
figure('Name', 'Análise Espectral (FFT)', 'Color', 'w', 'Position', [100 100 1200 800]);

% FFT do Sinal Bruto
fft_bruto = abs(fft(gyrox)) / N;
mag_bruto = 2 * fft_bruto(half_idx);

% FFT do Movimento Voluntário (CDF)
fft_cdf = abs(fft(x_cdf)) / N;
mag_cdf = 2 * fft_cdf(half_idx);

% FFT do WFLC
fft_wflc = abs(fft(tremor_wflc)) / N;
mag_wflc = 2 * fft_wflc(half_idx);

% FFT do BMFLC
fft_bmflc = abs(fft(tremor_bmflc)) / N;
mag_bmflc = 2 * fft_bmflc(half_idx);

% FFT do KF
fft_kf = abs(fft(tremor_kf)) / N;
mag_kf = 2 * fft_kf(half_idx);



idx_tremor = f_axis >= 4 & f_axis <= 15;
plot(f_axis(idx_tremor), mag_residual(idx_tremor), 'k--', 'LineWidth', 1.5); hold on;
plot(f_axis(idx_tremor), mag_wflc(idx_tremor), 'r', 'LineWidth', 1.5);
plot(f_axis(idx_tremor), mag_bmflc(idx_tremor), 'b', 'LineWidth', 1.5);
plot(f_axis(idx_tremor), mag_kf(idx_tremor), 'g', 'LineWidth', 1.5);
title('FFT da faixa de frequência do Tremor (4-15 Hz)');
xlabel('Frequência (Hz)'); ylabel('Magnitude');
xlim([4 15]); grid on;
legend('Resíduo', 'WFLC', 'BMFLC', 'KF', 'Location', 'best');

% Demonstração do Efeito da Compensação de Atraso (FMSEd)
figure('Name', 'Demonstração da Compensação de Atraso - FMSEd', 'Color', 'w', 'Position', [350 350 1200 800]);

% Para WFLC 
subplot(3,1,1);
plot(t, tremor_entrada, 'k', 'LineWidth', 1.5); hold on;
plot(t, tremor_wflc, 'r', 'LineWidth', 1);
title('WFLC: Sinal Original vs Estimado');
ylabel('Amplitude');
legend('Referência', 'Estimado');
grid on; axis tight; % Ajusta o eixo automaticamente ao tempo total

% Para BMFLC
subplot(3,1,2);
plot(t, tremor_entrada, 'k', 'LineWidth', 1.5); hold on;
plot(t, tremor_bmflc, 'b', 'LineWidth', 1);
title('BMFLC: Sinal Original vs Estimado');
ylabel('Amplitude');
legend('Referência', 'Estimado');
grid on; axis tight;

% Para KF
subplot(3,1,3);
plot(t, tremor_entrada, 'k', 'LineWidth', 1.5); hold on;
plot(t, tremor_kf, 'g', 'LineWidth', 1);
title('WFLC em cascata com o Filtro de Kalman: Sinal Original vs Estimado');
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Referência', 'Estimado');
grid on; axis tight;

%% FILTRO DE KALMAN: FUSÃO SENSORIAL PARA POSIÇÃO LINEAR

ax_bruto = dados(:, 2); 

ax_bruto = ax_bruto * 9.81; % Converte de 'g' para m/s² se necessário

% Filtro Passa-Banda no Acelerômetro (banda do tremor 4-15 Hz)
ax_tremor = bandpass(ax_bruto, [4 15], Fs);

% Estimativa do Raio (r = a / alpha)
% Derivada da velocidade angular (tremor_kf) para obter aceleração angular
alpha_tremor = diff(tremor_kf) ./ Ts;
alpha_tremor = [alpha_tremor(1); alpha_tremor]; 

% Cálculo automático do raio usando instantes de sinal forte
idx_v = abs(alpha_tremor) > max(abs(alpha_tremor)) * 0.2;
raio_est = mean(abs(ax_tremor(idx_v)) ./ abs(alpha_tremor(idx_v)));
raio_final = max(0.05, min(0.50, raio_est)); % Limita entre 5cm e 50cm

fprintf('--- Parâmetros de Entrada ---\n');
fprintf('Frequência: %.2f Hz | Raio Estimado: %.2f cm\n', Fs, raio_final*100);

%% FILTRO DE KALMAN DE FUSÃO
% Estados: x(1) = posição (m), x(2) = velocidade (m/s)

% Inicialização
x_f = [0; 0];
P_f = eye(2);
pos_linear_m = zeros(length(t), 1);

% Matrizes de Espaço de Estados
A = [1 Ts; 0 1];           % Transição: P = P0 + V*Ts
B = [0.5*Ts^2; Ts];        % Entrada: Efeito da aceleração
H = [0 1];                 % Medição: O giroscópio nos dá a Velocidade (v = w*r)

% Matrizes de Covariância (Ajuste de Confiança)
Q = [1e-7 0; 0 1e-5];      % Ruído de processo (confiança no modelo/acc)
R = 0.02;                  % Ruído de medição (confiança na ref do gyro)

for k = 1:length(t)
    % Entradas
    u = ax_tremor(k);                  % Aceleração (m/s²)
    v_ref = tremor_kf(k) * raio_final; % Velocidade de referência via Gyro (m/s)

    % --- PREDICÃO ---
    x_pred = A * x_f + B * u;
    P_pred = A * P_f * A' + Q;

    % --- CORREÇÃO ---
    S = H * P_pred * H' + R;
    K = P_pred * H' / S; % Ganho de Kalman

    x_f = x_pred + K * (v_ref - H * x_pred);
    P_f = (eye(2) - K * H) * P_pred;

    % Armazenar posição
    pos_linear_m(k) = x_f(1);
end



%% TRATAMENTO FINAL E CONVERSÃO
% Remover drift acumulado e converter metros para milímetros (mm)
pos_linear_mm = highpass(pos_linear_m, 0.5, Fs) * 1000;

%% PLOTS DE RESULTADOS
figure('Name', 'Fusão Sensorial: Deslocamento Linear do Tremor', 'Color', 'w', 'Position', [100 100 900 600]);
plot(t, pos_linear_mm, 'k', 'LineWidth', 1.5);
title('Posição Linear Estimada');
ylabel('Deslocamento (mm)');
xlabel('Tempo (s)');
xlim ([0 10]);
grid on;

% Estatísticas Rápidas
fprintf('\n--- Estatísticas do Deslocamento ---\n');
fprintf('Amplitude Pico-a-Pico: %.2f mm\n', (max(pos_linear_mm) - min(pos_linear_mm)));
fprintf('RMS do Tremor: %.2f mm\n', rms(pos_linear_mm));

referencia_tremor = pos_linear_mm;
tempo_referencia = t;
%% SALVAR REFERÊNCIA PARA O ESP32


dados_ref = referencia_tremor;

%%writematrix(dados_ref, 'ref.txt', 'Delimiter', 'tab');

%%fprintf('Referência para ESP32 salva com sucesso!\n');