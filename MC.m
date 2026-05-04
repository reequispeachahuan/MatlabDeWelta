clear all;
close all;
clc;
%%% Parámetros %%%  
Nt = 2; % Número de antenas transmisoras
Nr = 2; % Número de antenas receptoras
alpha = 2 * Nr; % Tamaño deseado de la red de comparadores
SNR_dB = -30:10:40; % Valores de SNR en dB
SNR = 10.^(SNR_dB / 10); % SNR en escala lineal
sigma_x = 1; % Energía total transmitida
channel_realizations = 1; % Número de realizaciones del canal
% Número de comparadores
full = Nr*(2*Nr-1); % Tamaño deseado de la red de comparadores
M_prime_full = 2 * Nr + full;
M_prime_random = 2 * Nr + alpha;
% Identidades y Covarianza
I_Nr_r = eye(2 * Nr);
Cx_r = (1/2) * sigma_x^2 * eye(2 * Nt);
% Inicializar capacidades
capacities_optimized = zeros(length(SNR), channel_realizations);
capacities_montecarlo = zeros(length(SNR), channel_realizations); % Capacidades optimizadas con Monte Carlo
I_random = zeros(length(SNR), channel_realizations); % Capacidades de la red aleatoria
I_full = zeros(length(SNR), channel_realizations); % Capacidades de la red completa
I_withoutB = zeros(length(SNR), channel_realizations); % Capacidades de la red sin comparadores
%%%  
B_alpha_f = 1/sqrt(2) * create_comparator_matrix(Nr, 'full', []);
B_full = [I_Nr_r ; B_alpha_f];
% Generamos los comparadores optimizados
%var_H_r = sum(abs(B_prime).^2, 2); % Calculamos la varianza
%[~, position] = sort(var_H_r, 'descend'); % Ordenamos la varianza de H_r

% Inicialización de la barra de progreso
h = waitbar(0, 'Comenzando simulación...');
for i_channel = 1:channel_realizations
    for i = 1:length(SNR)
        %% Generación del canal H
        H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);
        H_r = [real(H), -imag(H); imag(H), real(H)];
        
        %% Cálculo de la covarianza
        sigma_n = sqrt(sigma_x^2 / SNR(i)); 
        Cn_r = (sigma_n^2 / 2) * I_Nr_r; % Covarianza del ruido (real)
        %% Matriz de covarianza
        Cz_r_full = B_full * (H_r * Cx_r * H_r') * B_full' + B_full * Cn_r * B_full';
        
        %% Normalización y cálculo de H_eff_r_q
        lambda = (2 / pi) * ((pi / 2 - 1) + (sigma_n^2 / (2 * (Nt * sigma_x^2 / 2 + sigma_n^2 / 2))));
        k_r_full = diag(1 ./ sqrt(diag(B_full * (H_r * Cx_r * H_r') * B_full' + B_full * (sigma_n^2 / 2) * I_Nr_r * B_full')));
        H_eff_r_q_full = sqrt(2 / pi) * k_r_full * B_full * H_r;
        
        [Deltao, ~] = solve_selection_yalmip(M_prime_full, Nr, alpha, Nt, lambda, sigma_x, H_eff_r_q_full);
        
        % Seleccionar los índices de los comparadores
        [~, sorted_indices] = maxk(Deltao, 2 * Nr + alpha);
        vector_delta_0 = zeros(M_prime_full, 1);
        vector_delta_0(sorted_indices) = 1;
        selected_indices = find(vector_delta_0);  
        B_select_alpha = B_full(selected_indices, :);
        
        %% Cálculos de capacidad
        Cz_r = B_select_alpha * (H_r * Cx_r * H_r') * B_select_alpha' + B_select_alpha * Cn_r * B_select_alpha';
        k_r = diag(1 ./ sqrt(diag(Cz_r)));
        H_eff_r_q1 = sqrt(2 / pi) * k_r * B_select_alpha * H_r;
        
        C_eta_eff_r = (2 / pi) * (asin(k_r * Cz_r * k_r) - k_r * Cz_r * k_r) + k_r * B_select_alpha * Cn_r * B_select_alpha' * k_r;
        capacities_optimized(i, i_channel) = 1 / 2 * log2(det(eye(2 * Nr + alpha) + pinv(real(C_eta_eff_r)) * ((sigma_x^2 / 2) * (H_eff_r_q1 * H_eff_r_q1'))));
        %% Monte Carlo - Selección Aleatoria de Comparadores
        num_simulations = 1; % Número de simulaciones Monte Carlo
        best_capacity_montecarlo = -inf; % Inicializar la mejor capacidad para Monte Carlo
        for sim = 1:num_simulations
            % Generación aleatoria de una configuración de comparadores
            random_indices = randperm(M_prime_full, 2 * Nr + alpha); % Selección aleatoria de indices
            B_select_montecarlo = B_full(random_indices, :); % Red seleccionada aleatoriamente
            % Calcular la capacidad para la configuración actual
            Cz_r = B_select_montecarlo * (H_r * Cx_r * H_r') * B_select_montecarlo' + B_select_montecarlo * Cn_r * B_select_montecarlo';
            k_r = diag(1 ./ sqrt(diag(Cz_r)));
            H_eff_r_q2 = sqrt(2 / pi) * k_r * B_select_montecarlo * H_r;
            C_eta_eff_r = (2 / pi) * (asin(k_r * Cz_r * k_r) - k_r * Cz_r * k_r) + k_r * B_select_montecarlo * Cn_r * B_select_montecarlo' * k_r;
            capacity_montecarlo = 1 / 2 * log2(det(eye(2 * Nr + alpha) + pinv(real(C_eta_eff_r)) * ((sigma_x^2 / 2) * (H_eff_r_q2 * H_eff_r_q2'))));
            % Actualizar la mejor capacidad
            if capacity_montecarlo > best_capacity_montecarlo
                best_capacity_montecarlo = capacity_montecarlo;
            end
        end
        % Almacenar la capacidad obtenida con Monte Carlo
        capacities_montecarlo(i, i_channel) = best_capacity_montecarlo;
        %% Red Aleatoria - Capacidad
        B_random_indices = randperm(M_prime_random, 2 * Nr + alpha);
        B_random_select = B_full(B_random_indices, :); % Red aleatoria seleccionada
        Cz_r_random = B_random_select * (H_r * Cx_r * H_r') * B_random_select' + B_random_select * Cn_r * B_random_select';
        k_r_random = diag(1 ./ sqrt(diag(Cz_r_random)));
        H_eff_r_q_random = sqrt(2 / pi) * k_r_random * B_random_select * H_r; 
        C_eta_eff_r_random = (2 / pi) * (asin(k_r_random * Cz_r_random * k_r_random) - k_r_random * Cz_r_random * k_r_random) + k_r_random * B_random_select * Cn_r * B_random_select' * k_r_random;
        I_random(i, i_channel) = 1 / 2 * log2(det(eye(M_prime_random) + pinv(real(C_eta_eff_r_random)) * ((sigma_x^2 / 2) * (H_eff_r_q_random * H_eff_r_q_random'))));
        %% Red Completa - Capacidad
        Cz_r_full = B_full * (H_r * Cx_r * H_r') * B_full' + B_full * Cn_r * B_full';
        k_r_full = diag(1 ./ sqrt(diag(Cz_r_full)));
        H_eff_r_q_full = sqrt(2 / pi) * k_r_full * B_full * H_r; 
        C_eta_eff_r_full = (2 / pi) * (asin(k_r_full * Cz_r_full * k_r_full) - k_r_full * Cz_r_full * k_r_full) + k_r_full * B_full * Cn_r * B_full' * k_r_full;
        I_full(i, i_channel) = 1 / 2 * log2(det(eye(M_prime_full) + pinv(real(C_eta_eff_r_full)) * ((sigma_x^2 / 2) * (H_eff_r_q_full * H_eff_r_q_full'))));
        %% Red Sin Comparadores - Capacidad
        Cz_r_withoutB = (H_r * Cx_r * H_r') + (sigma_n^2 / 2) * I_Nr_r;
        k_r_withoutB = diag(1 ./ sqrt(diag(Cz_r_withoutB)));
        H_eff_r_q_withoutB = sqrt(2 / pi) * k_r_withoutB * H_r; 
        C_eta_eff_r_withoutB = (2 / pi) * (asin(k_r_withoutB * Cz_r_withoutB * k_r_withoutB) - k_r_withoutB * Cz_r_withoutB * k_r_withoutB) + k_r_withoutB * Cn_r * k_r_withoutB;
        I_withoutB(i, i_channel) = 1 / 2 * log2(det(eye(2 * Nr) + pinv(real(C_eta_eff_r_withoutB)) * ((sigma_x^2 / 2) * (H_eff_r_q_withoutB * H_eff_r_q_withoutB'))));
    end
    % Actualizar la barra de progreso
    waitbar(i_channel / channel_realizations, h, sprintf('Simulaciones: %d/%d', i_channel, channel_realizations));
end
% Cerrar la barra de progreso
close(h);
% Promediar las capacidades
capacities_optimized_av = mean(capacities_optimized, 2);
capacities_montecarlo_av = mean(capacities_montecarlo, 2);
I_random_av = mean(I_random, 2); % Promedio de la red aleatoria
I_full_av = mean(I_full, 2); % Promedio de la red completa
I_withoutB_av = mean(I_withoutB, 2); % Promedio de la red sin comparadores
%%% Graficar la capacidad optimizada para todas las redes %%% 
figure;
plot(SNR_dB, capacities_optimized_av, 'b-o', 'LineWidth', 2, 'DisplayName', 'Proposed Optimized Comparator Network');
hold on;
plot(SNR_dB, capacities_montecarlo_av, 'm-*', 'LineWidth', 2, 'DisplayName', 'Monte Carlo Optimized');
plot(SNR_dB, I_random_av, 'g-s', 'LineWidth', 2, 'DisplayName', 'Random Comparator Network');
plot(SNR_dB, I_full_av, 'r-x', 'LineWidth', 2, 'DisplayName', 'Full Comparator Network');
plot(SNR_dB, I_withoutB_av, 'c-d', 'LineWidth', 2, 'DisplayName', 'Without Comparator Network');
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Capacidad (bits/s/Hz)', 'Interpreter', 'latex');
grid on;
title('Capacidad optimizada vs. SNR', 'Interpreter', 'latex');
legend('show');
