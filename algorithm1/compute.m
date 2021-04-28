function [E, values] = compute(I,J,U,V)
X = U*V';
[row col] = size(X);
E = -X;
for k = 1:length(I)
    E(I(k),J(k)) = 0;
end
