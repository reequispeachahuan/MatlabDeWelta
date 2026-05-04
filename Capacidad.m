clear all;
close all;
clc;

%%% Parámetros %%%  
Nt = 2; % Número de antenas transmisoras
Nr = 2; % Número de antenas receptor
alpha = 2 * Nr; % Tamaño deseado de la red de comparadores (alpha = 8 * Nr)
SNR_dB = -30:10:40; % Valores de SNR en dB
SNR = 10.^(SNR_dB / 10); % SNR 
sigma_x = 1; % Energía total transmitida
channel_realizations = 1; % Número de realizaciones del canal
full = Nr*(2*Nr-1); % Tamaño deseado de la red de comparadores
M_prime_full = 2 * Nr + full;
M_prime_random = 2 * Nr + alpha;

% Identidades y Covarianza
I_Nr_r = eye(2 * Nr);
Cx_r = (1/2) * sigma_x^2 * eye(2 * Nt);

% Inicializar capacidades
capacities_optimized = zeros(length(SNR), 1);
capacities_random = zeros(length(SNR), 1);
I_full = zeros(length(SNR), 1); % Inicialización de la variable I_full

% Comparator Network%%  
B_prime = 1/sqrt(2) * create_comparator_matrix(Nr, 'partial', alpha);
B_random = [I_Nr_r ; B_prime];

% Fully
B_alpha_f = 1/sqrt(2) * create_comparator_matrix(Nr, 'full', []);
B_full = [I_Nr_r ; B_alpha_f];

% Inicialización de la capacidad óptima (búsqueda exhaustiva)
capacities_exhaustive = zeros(length(SNR), 1);

for i_channel = 1:channel_realizations
    for i = 1:length(SNR)
        %% Generación del canal H
        H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);
        H_r = [real(H), -imag(H); imag(H), real(H)];

        %% Cálculo de la covarianza
        sigma_n = sqrt(sigma_x^2 / SNR(i)); 
        Cn_r = (sigma_n^2 / 2) * I_Nr_r; % Covarianza del ruido (real)

        % Matriz de covarianza cuantizada exacta
        Cz_r_full = B_full * (H_r * Cx_r * H_r') * B_full' + B_full * Cn_r * B_full';

        lambda = (2 / pi)*((pi/2 - 1) + (sigma_n^2 / (2 * (Nt * sigma_x^2/2 + sigma_n^2 / 2)))); 
        k_r_full = diag(1 ./ sqrt(diag(B_full * (H_r * Cx_r * H_r') * B_full' + B_full * (sigma_n^2 / 2) * I_Nr_r * B_full'))); 
        H_eff_r_q_full = sqrt(2 / pi) * k_r_full * B_full * H_r;

        % Aplica la transposición correctamente para H_eff_r_q
        H_eff_r_q_transpose = H_eff_r_q_full.'; 

        % Optimización Convexa con CVX para Selección Completa
        cvx_begin quiet sdp
            variable Delta(M_prime_full); 
            % Maximizar capacidad
            maximize(1/2 * log_det( eye(2*Nt) + 1/lambda *(sigma_x^2/2) * ((H_eff_r_q_full' * diag(Delta) * H_eff_r_q_full)) ));
            % Restricciones
            subject to
                for i_delta = 1:2*Nr
                    Delta(i_delta) == 1;
                end
                0 <= Delta <= 1;
                sum(Delta) == 2*Nr + alpha;
        cvx_end

        % Seleccionar los índices de los comparadores
        [~, sorted_indices] = maxk(Delta, 2 * Nr + alpha);
        vector_delta_0 = zeros(M_prime_full, 1);
        vector_delta_0(sorted_indices) = 1;
        selected_indices = find(vector_delta_0);  
        B_select_alpha = B_full(selected_indices, :);

        % Cálculos de capacidad para la red optimizada
        Cz_r = B_select_alpha * (H_r * Cx_r * H_r') * B_select_alpha' + B_select_alpha * Cn_r * B_select_alpha';
        k_r = diag(1 ./ sqrt(diag(Cz_r))); 
        H_eff_r_q1 = sqrt(2 / pi) * k_r * B_select_alpha * H_r;

        C_eta_eff_r = (2 / pi) * (asin(k_r * Cz_r * k_r) - k_r * Cz_r * k_r) + k_r * B_select_alpha * Cn_r * B_select_alpha' * k_r;
        I_select(i, i_channel) = 1/2 * log2(det(eye(2*Nr+alpha) + pinv(real(C_eta_eff_r)) * (((sigma_x^2/2) * (H_eff_r_q1*H_eff_r_q1')))));
        %
        % Please compute here H_eff_r_q and C_eta`R   with B_full
     C_eta_eff_r_full = (2 / pi) * (asin(k_r_full * Cz_r_full * k_r_full) - k_r_full * Cz_r_full * k_r_full) + k_r_full * B_full * Cn_r * B_full' * k_r_full;
% 
I_full(i,i_channel) = 1/2 * log2(det(eye(M_prime_full) + pinv(real(C_eta_eff_r_full)) * (((sigma_x^2/2)* (H_eff_r_q_full*H_eff_r_q_full')))));
% 
        % Red aleatoria
       % B_prime = 1/sqrt(2) * get_random_perm(alpha, 2 * Nr);
        %B_random = [I_Nr_r; B_prime];
        Cz_r_random = B_random * (H_r * Cx_r * H_r') * B_random' + B_random * (sigma_n^2 / 2) * I_Nr_r * B_random';
        k_r_random = diag(1 ./ sqrt(diag(Cz_r_random))); 
        H_eff_r_q_random = sqrt(2 / pi) * k_r_random * B_random * H_r; 
        C_eta_eff_r_random = (2 / pi) * (asin(k_r_random * Cz_r_random * k_r_random) - k_r_random * Cz_r_random * k_r_random) + k_r_random * B_random * Cn_r * B_random' * k_r_random;
        I_random(i, i_channel) = 1/2 * log2(det(eye(M_prime_random) + pinv(real(C_eta_eff_r_random)) * (((sigma_x^2/2)* (H_eff_r_q_random*H_eff_r_q_random')))));

        % Red sin comparadores
        Cz_r = (H_r * Cx_r * H_r') + (sigma_n^2 / 2) * I_Nr_r;
        k_r = diag(1 ./ sqrt(diag(Cz_r))); 
        H_eff_r_q = sqrt(2 / pi) * k_r * H_r; 
        C_eta_eff_r = (2 / pi) * (asin(k_r * Cz_r * k_r) - k_r * Cz_r * k_r) + k_r * Cn_r * k_r;
        I_withoutB(i, i_channel) = 1/2 * log2(det(eye(2*Nr) + pinv(real(C_eta_eff_r)) * (((sigma_x^2/2)* (H_eff_r_q*H_eff_r_q')))));

        % Búsqueda Exhaustiva: Evaluar todas las combinaciones posibles
        k = 2 * Nr + alpha;  % Número de comparadores seleccionados
        combinations = nchoosek(1:M_prime_full, k);  % Obtener todas las combinaciones posibles
        max_capacity = -inf;  % Inicializamos la capacidad máxima

        % Evaluamos todas las combinaciones
        for j = 1:size(combinations, 1)
            selected_indices_exh = combinations(j, :);
            B_select_exh = B_full(selected_indices_exh, :);
            Cz_r_exh = B_select_exh * (H_r * Cx_r * H_r') * B_select_exh' + B_select_exh * Cn_r * B_select_exh';
            k_r_exh = diag(1 ./ sqrt(diag(Cz_r_exh))); 
            H_eff_r_q_exh = sqrt(2 / pi) * k_r_exh * B_select_exh * H_r;

            C_eta_eff_r_exh = (2 / pi) * (asin(k_r_exh * Cz_r_exh * k_r_exh) - k_r_exh * Cz_r_exh * k_r_exh) + k_r_exh * B_select_exh * Cn_r * B_select_exh' * k_r_exh;
            capacity = 1/2 * log2(det(eye(2*Nr+alpha) + pinv(real(C_eta_eff_r_exh)) * (((sigma_x^2/2) * (H_eff_r_q_exh*H_eff_r_q_exh')))));

            % Actualizar la capacidad máxima
            if capacity > max_capacity
                max_capacity = capacity;
            end
        end
        capacities_exhaustive(i, i_channel) = max_capacity;
    end
end

% Promediar las capacidades
I_random_av = sum(I_random, 2) / channel_realizations;
I_select_av = sum(I_select, 2) / channel_realizations;
I_full_av = sum(I_full, 2) / channel_realizations; % Promedio de la capacidad de la red completa
I_withoutB_av = sum(I_withoutB, 2) / channel_realizations;
I_exhaustive_av = sum(capacities_exhaustive, 2) / channel_realizations;

%%% %%% 
figure;
plot(SNR_dB, I_random_av, 'g-s', 'LineWidth', 2, 'DisplayName', 'Random Comparator Network');
hold on;
plot(SNR_dB, I_select_av, 'b-o', 'LineWidth', 2, 'DisplayName', 'Proposed Optimized Comparator Network');
plot(SNR_dB, I_full_av, 'r-x', 'LineWidth', 2, 'DisplayName', 'Full Comparator Network');
plot(SNR_dB, I_withoutB_av, 'c-x', 'LineWidth', 2, 'DisplayName', 'Without Comparator Network');

% 
plot(SNR_dB, I_exhaustive_av, 'k-^', 'LineWidth', 2, 'DisplayName', 'Optimal (Exhaustive Search)'); 

xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Capacidad (bits/s/Hz)', 'Interpreter', 'latex');
grid on;
title('Capacidad optimizada vs. SNR', 'Interpreter', 'latex');
legend('show');