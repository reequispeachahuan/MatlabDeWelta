function m_comb = get_all_perm(n_comb, cols)

m_comb = zeros([n_comb 2]);

row = 0;
for i=1:cols
    for j=i+1:cols
        m_comb(row+1, 1) = i;
        m_comb(row+1, 2) = j;
        row = row + 1;
    end
end
end
