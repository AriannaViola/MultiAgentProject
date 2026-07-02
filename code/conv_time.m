function tcv = conv_time(d, t, tol)
%CONV_TIME  Primo istante in cui il disaccordo scende sotto tol (NaN se mai).
idx = find(d < tol, 1, 'first');
if isempty(idx), tcv = NaN; else, tcv = t(idx); end
end
