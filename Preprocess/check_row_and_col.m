function satisfy = check_row_and_col(rows,cols,selected_rows,selected_cols)
    row = zeros(1,rows);
    col = zeros(1,cols);
    for i = 1:rows
        row(i) = ismember(i,selected_rows);
    end
    for i = 1:cols
        col(i) = ismember(i,selected_cols);
    end
    satisfy = (all(row)& all(col));
end