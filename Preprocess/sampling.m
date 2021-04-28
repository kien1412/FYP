function [observed_entries, observed_values] = sampling(x)
    [rows, cols] = size(x);
    number_of_entries = rows*cols;
    %sampling_matrix = zeros(rows, cols);
    sample_size = floor(0.25*number_of_entries); %sample size as percentage of total entries
    observed_values = zeros(1,sample_size);
    
    observed_entries = randsample(1:number_of_entries, sample_size);
    [selected_row, selected_col] = ind2sub([rows cols], observed_entries);
    for i = 1:sample_size
        observed_values(i) = x(selected_row(i), selected_col(i));
    end
    
    while ~check_row_and_col(rows,cols,selected_row,selected_col)
        observed_entries = randsample(1:number_of_entries, sample_size);
        [selected_row, selected_col] = ind2sub([rows cols], observed_entries);
        for i = 1:sample_size
            observed_values(i) = x(selected_row(i), selected_col(i));
        end
    end
end

