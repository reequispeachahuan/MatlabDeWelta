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
channel_realizations = 20; % Número de realizaciones del canal
full = Nr*(2*Nr-1); % Tamaño deseado de la red de comparadores
% Número de comparadores
M_prime_full = 2 * Nr + full;
M_prime_random = 2 * Nr + alpha;
% Identidades y Covarianza
I_Nr_r = eye(2 * Nr);
Cx_r = (1/2) * sigma_x^2 * eye(2 * Nt);
% Inicializar capacidades
capacities_optimized = zeros(length(SNR), 1);
capacities_random = zeros(length(SNR), 1);
%%% Formación de la Red de Comparadores %%%  
B_prime = 1/sqrt(2) * get_random_perm(alpha, 2 * Nr); % Comparadores aleatorios generados aleatoriamente
B_random = [I_Nr_r; B_prime];
% Generamos los comparadores optimizados
var_H_r = sum(abs(B_prime).^2, 2); % Calculamos la varianza
[~, position] = sort(var_H_r, 'descend'); % Ordenamos la varianza de H_r
% Generación de los comparadores optimizados
B_alpha_f = 1/sqrt(2) * get_random_perm(full, 2 * Nr); % Comparadores optimizados basados en la permutación
B_full = [I_Nr_r; B_alpha_f]; % Formar B_full al añadir I_Nr_r y los comparadores optimizados
for i_channel=1:channel_realizations
%%% Optimización Convexa con CVX para Selección Completa %%%  
for i = 1:length(SNR)
    %% 
    H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);
    H_r = [real(H), -imag(H); imag(H), real(H)];
    %%
    sigma_n = sqrt(sigma_x^2 / SNR(i)); 
    Cn_r = (sigma_n^2 / 2) * I_Nr_r; % Covarianza del ruido (real)
    % Matriz de covarianza cuantizada exacta
    Cz_r_full = B_full * (H_r * Cx_r * H_r') * B_full' + B_full * Cn_r * B_full';
    % 
    lambda = (2 / pi)*((pi/2 - 1) + (sigma_n^2 / (2 * (Nt * sigma_x^2/2 + sigma_n^2 / 2)))); 
    % Normalización y cálculo de H_eff_r_q
    k_r_full = diag(1 ./ sqrt(diag(B_full * (H_r * Cx_r * H_r') * B_full' + B_full * (sigma_n^2 / 2) * I_Nr_r * B_full'))); 
    H_eff_r_q_full= sqrt(2 / pi) * k_r_full * B_full * H_r;
    % Aplica la transposición correctamente para H_eff_r_q
    H_eff_r_q_transpose = H_eff_r_q_full.'; % Transposición de H_eff_r_q
    I_matrix = eye(2 * Nt);
    
    cvx_begin quiet sdp
        variable Deltao(M_prime_full); 
        % Maximizar capacidad
        maximize(1/2 * log_det( eye(2*Nt) + 1/lambda *(sigma_x^2/2) * ((H_eff_r_q_full' * diag(Deltao) * H_eff_r_q_full))));
        % Restricciones
        subject to
            for i_delta = 1:2*Nr
                Deltao(i_delta) == 1;
            end
            0 <= Deltao <= 1;
            sum(Deltao) == 2*Nr + alpha;
    cvx_end
    % Seleccionar los índices de los comparadores
    [~, sorted_indices] = maxk(Deltao, 2 * Nr + alpha);
    vector_delta_0 = zeros(M_prime_full, 1);
    vector_delta_0(sorted_indices) = 1;
    % 
    selected_indices = find(vector_delta_0);  % 
B_select_alpha = B_full(selected_indices   , :   );
       Cz_r = B_select_alpha * (H_r * Cx_r * H_r') * B_select_alpha' + B_select_alpha * Cn_r * B_select_alpha';
        k_r = diag(1 ./ sqrt(diag(Cz_r))); % 
        H_eff_r_q1 = sqrt(2 / pi) * k_r * B_select_alpha * H_r;
        %
        C_eta_eff_r = (2 / pi) * (asin(k_r * Cz_r * k_r) - k_r * Cz_r * k_r) + k_r * B_select_alpha * Cn_r * B_select_alpha' * k_r;
% Please compute here the true capacity (14) with the selected comparator
% netwoek
I_select(i,i_channel) = 1/2 * log2(det(eye(2*Nr+alpha) + pinv(real(C_eta_eff_r)) * (((sigma_x^2/2)* (H_eff_r_q1*H_eff_r_q1')))));
% Please compute here H_eff_r_q and C_eta`R   with B_full
C_eta_eff_r_full = (2 / pi) * (asin(k_r_full * Cz_r_full * k_r_full) - k_r_full * Cz_r_full * k_r_full) + k_r_full * B_full * Cn_r * B_full' * k_r_full;
% Please compute here the true capacity (14) with the selected comparator
% netwoek
I_full(i,i_channel) = 1/2 * log2(det(eye(M_prime_full) + pinv(real(C_eta_eff_r_full)) * (((sigma_x^2/2)* (H_eff_r_q_full*H_eff_r_q_full')))));
% Please compute here the B_alpha with alpha random comparators
%B_alpha== 1/sqrt(2) * get_random_perm(alpha, 2 * Nr);
% Please compute here the B_random
%%% Formación de la Red de Comparadores %%%  
B_prime = 1/sqrt(2) * get_random_perm(alpha, 2 * Nr); % Comparadores aleatorios generados aleatoriamente
B_random = [I_Nr_r; B_prime];
%%%
        Cz_r_random = B_random * (H_r * Cx_r * H_r') * B_random' + B_random * (sigma_n^2 / 2) * I_Nr_r * B_random';
        k_r_random = diag(1 ./ sqrt(diag(Cz_r_random))); 
        H_eff_r_q_random = sqrt(2 / pi) * k_r_random * B_random * H_r; 
        C_eta_eff_r_random = (2 / pi) * (asin(k_r_random * Cz_r_random * k_r_random) - k_r_random * Cz_r_random * k_r_random) + k_r_random * B_random * Cn_r * B_random' * k_r_random
        I_random(i,i_channel) = 1/2 * log2(det(eye(M_prime_random) + pinv(real(C_eta_eff_r_random)) * (((sigma_x^2/2)* (H_eff_r_q_random*H_eff_r_q_random')))));
% Please compute here H_eff_r_q and C_eta`R   without B
        Cz_r= (H_r * Cx_r * H_r')  +  (sigma_n^2 / 2) * I_Nr_r;
        k_r= diag(1 ./ sqrt(diag(Cz_r))); 
        H_eff_r_q= sqrt(2 / pi) * k_r * H_r; 
        C_eta_eff_r = (2 / pi) * (asin(k_r* Cz_r * k_r) - k_r * Cz_r * k_r) + k_r * Cn_r* k_r
        I_withoutB(i,i_channel) = 1/2 * log2(det(eye(2*Nr) + pinv(real(C_eta_eff_r)) * (((sigma_x^2/2)* (H_eff_r_q*H_eff_r_q')))));
    % Seleccionar las filas de H
    %H_selected = H(selected_indices, :);  % Selecciona las filas correspondientes a los índices seleccionados
    % Calcular la capacidad con las filas seleccionadas
    %capacity_selected = log2(det(eye(length(selected_indices)) + (H_selected' * H_selected)));
    % Almacenar la capacidad optimizada
    capacities_optimized(i) = cvx_optval;
    
    % Mostrar resultados
    disp(['SNR = ', num2str(SNR_dB(i)), ' dB']);
    disp('Filas seleccionadas (Selección Completa):');
    disp(selected_indices);
    % Crear la nueva matriz B seleccionada solo con las filas seleccionadas
    B_selected = B_full(selected_indices, :); % 
    % Mostrar la submatriz B seleccionada
    disp('Submatriz B seleccionada (con las filas seleccionadas):');
    disp(B_selected);
end
end
 I_random_av  = sum(I_random, 2) / channel_realizations;
 I_select_av  = sum(I_select, 2) / channel_realizations;
 I_full_av  = sum(I_full, 2) / channel_realizations;
 I_withoutB_av=sum(I_withoutB, 2) / channel_realizations;
%%% Graficar la capacidad optimizada para ambas selecciones %%% 
%figure; 
%plot(SNR_dB, capacities_optimized, 'r-s', 'LineWidth', 2, 'DisplayName', 'Capacidad optimizada (Full Selection)'); 
%hold on; 
%ylabel('Capacidad (bits/s/Hz)');
%grid on;
%title('Capacidad optimizada vs. SNR para Selección Completa y Aleatoria', 'Interpreter', 'latex');
%legend('show');
%%% 
figure;
plot(SNR_dB, I_random_av, 'g-s', 'LineWidth', 2, 'DisplayName', 'Random Comparator Network');
hold on;
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Capacidad (bits/s/Hz)', 'Interpreter', 'latex');
plot(SNR_dB, I_select_av, 'b-o', 'LineWidth', 2, 'DisplayName', 'Proposed Optimized Comparator Network');
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Capacidad (bits/s/Hz)', 'Interpreter', 'latex');
plot(SNR_dB, I_full_av, 'r-x', 'LineWidth', 2, 'DisplayName', 'Full Comparator Network');
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Capacidad (bits/s/Hz)', 'Interpreter', 'latex')
plot(SNR_dB, I_withoutB_av, 'c-x', 'LineWidth', 2, 'DisplayName', 'Without Comparator Network');
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Capacidad (bits/s/Hz)', 'Interpreter', 'latex')
grid on;
title('Capacidad optimizada vs. SNR', 'Interpreter', 'latex');
legend('show');