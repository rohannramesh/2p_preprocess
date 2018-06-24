function var = remove_nans(var)

if sum(size(var) == 1)
    var(isnan(var)) = [];
else
    ind_nan = find(isnan(var(:,1)));
    var(ind_nan,:) = [];
end