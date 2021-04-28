function noisy_matrix = add_noise(x)
    [rows, cols] = size(x);
    number_of_entries = rows*cols;
    %variance = sum(x, 'all')/number_of_entries;
    variance = 0.01;
    noisy_sample = 0.50*number_of_entries;
    noisy_entries = randsample(1:number_of_entries, noisy_sample);
    [selected_row, selected_col] = ind2sub([rows cols], noisy_entries);
    noisy_matrix = x;
    for i = 1:noisy_sample
        noisy_matrix(selected_row(i), selected_col(i)) = noisy_matrix(selected_row(i), selected_col(i)) + (variance.^0.5)*randn();
    end
end
    