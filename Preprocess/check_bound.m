function [ratio, lower_bound, upper_bound] = check_bound(X, k)
    [rows,cols] = size(X);
    %observed_entries = round(0.4*k*(rows+cols)*log(rows+cols));
    n = max([rows,cols]);
    total_entries = rows*cols;
    observed_entries = round(1.1*(4*n*k - 4*k*k));
    ratio = observed_entries/total_entries;
    lower = 4*n*k - 4*k*k;
    upper = rows*cols;
    lower_bound = (observed_entries >= lower);
    upper_bound = (observed_entries < upper);
end
    
    