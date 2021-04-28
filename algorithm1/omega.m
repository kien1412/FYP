function values = omega(I,J,U,V)
p = length(I);
X = U*V';
values = zeros(p, 1);
for k = 1:p
    values(k) = X(I(k),J(k));
end
