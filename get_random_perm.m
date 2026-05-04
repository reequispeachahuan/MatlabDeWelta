function matrix_final = get_random_perm(n_rows, n_cols)

% random number generator
rng('shuffle');

% max number of possible combinations for matrix B
n_max_comb = 0;
for i=1:n_cols-1
    n_max_comb = n_max_comb + i;
end

% [ [1,2], [1,3], ..., [1, col], ..., [2,3], [2,4], ..., [2, col], ...
% [col-1, col]
matrix_all_perm_indexes = get_all_perm(n_max_comb, n_cols);

% choosing n_rows random indexes from 1 to n_max_comb
array_random_rows_indexes = randperm(n_max_comb, n_rows);

% select just the combinations of the lines sorted above
m_with_selected_perm_indexes = matrix_all_perm_indexes(array_random_rows_indexes, :);

% finding the resulting matrix
matrix_final = zeros([n_rows, n_cols]);
for i=1:n_rows
    % setting 1 in the selected indexes
    matrix_final(i,m_with_selected_perm_indexes(i,1)) = 1;
    matrix_final(i,m_with_selected_perm_indexes(i,2)) = 1;
    
    % randomly choosing which one wil be -1
    negative_index = randi(2);
    matrix_final(i,m_with_selected_perm_indexes(i,negative_index)) = -1;
end

end