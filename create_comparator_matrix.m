function B_prime = create_comparator_matrix(Nr, type, num_comparators)
    % Nr: Número de antenas receptoras.
    % type: 'full' para red totalmente conectada, 'partial' para red parcialmente conectada.
    % num_comparators: Número de comparadores deseados (solo aplica si type='partial').

    % Número total de columnas (2Nr debido a la extensión real-imaginaria)
    n_cols = 2 * Nr;

    % Número total de combinaciones posibles para red totalmente conectada
    alpha_full = Nr * (2 * Nr - 1);

    % Validar el tipo de red
    if strcmp(type, 'full')
        % Usar todas las combinaciones posibles
        combinations = nchoosek(1:n_cols, 2);
        selected_combinations = combinations;
    elseif strcmp(type, 'partial')
        % Validar que num_comparators no exceda el máximo permitido
        if num_comparators > alpha_full
            error('El número de comparadores excede el máximo permitido.');
        end
        combinations = nchoosek(1:n_cols, 2);
        selected_indices = randperm(size(combinations, 1), num_comparators);
        selected_combinations = combinations(selected_indices, :);
    else
        error('Tipo de red inválido. Usa "full" o "partial".');
    end

    % Crear matriz B_prime
    n_rows = size(selected_combinations, 1); % Número de filas según las combinaciones seleccionadas
    B_prime = zeros(n_rows, n_cols); % Inicializar matriz

    for i = 1:n_rows
        %
        B_prime(i, selected_combinations(i, 1)) = 1;
        % 
        B_prime(i, selected_combinations(i, 2)) = -1;
    end
end