function largest_value = largest_value(matrix)
    [~,ii] = max(abs(matrix));
    largest_value = matrix(ii);
end