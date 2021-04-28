function [A, primal_residual, muk, error, rak] = algorithm_2(observed_entries, observed_values, m, n, tol, maxIter, X_org, noise_tol)
global A Sparse_Z;
%addpath PROPACK;


% read sparse matrix D
%[m,n] = size(D);
%[I,J,V] = find(D);
[I,J] = ind2sub([m n], observed_entries);
p = length(I);
I = reshape(I, [p,1]);
J = reshape(J, [p,1]);
V = observed_values;
V = reshape(V, [p,1]);

noise_tol = sqrt((p + sqrt(8*p))*noise_tol);
%col = [0; find(diff(J)); p];
D = zeros(m,n);
for i = 1:p
    D(I(i),J(i)) = V(i);
end
Sparse_Z = D;


% initialize all 0

E = zeros(m,n);
X = zeros(m,n);
Y_1 = sparse(m,n);
A.U = zeros(m, 5);
A.V = zeros(n, 5);
d_norm = norm(V, 'fro');
%mu = 1/(lansvd('Axz','Atxz',m,n,1,'L'));
mu = 1/svds(D,1);
rho_s = p/(m*n);
rho = 1.2172 + 1.8588*rho_s;
sv = 5;
svp = sv;


% Iteration 
iter = 0;
converged = false;
stopCriterion = 1;
%error_1 = zeros(1,500);
primal_residual = zeros(1,maxIter);
muk = zeros(1,maxIter+1);
error = zeros(1,maxIter);
rak = zeros(1,maxIter);

muk(iter + 1) = mu;
while ~converged         
    %% alternative projection 
    iter = iter + 1;    
    Sparse_Z = D - E - 1/mu*Y_1;
    
    if stopCriterion > 10*tol
        options.tol = 10*tol;
    else
        options.tol = min(0.1*tol, 0.01/mu);
    end
    %[A.U,S,A.V] = lansvd('Axz','Atxz',m,n,sv,'L',options); 
    %[q,w,e] = find(Sparse_Z);
    %pro(iter) = length(e)/(m*n);
    %rak_m(iter) = rank(Sparse_Z);
    
    [A.U,S,A.V] = svds(Sparse_Z, sv);
    %[A.U,S,A.V] = lansvd(Sparse_Z,sv,'L',options);
    %% predict the rank of A.
    diagS = diag(S);
    diagS = diagS(1:sv);
    svn = length(find(diagS > 1/mu));
    svp = svn;
    ratio = diagS(2:end-1)./diagS(3:end);
    [max_ratio, max_idx] = max(ratio);
    if max_ratio > 2
        svp = min(svn, max_idx + 1);
    end
    if svp < sv 
        sv = min(svp + 1, m);
    else
        sv = min(svp + 10, m);
    end
        
    %% update A Y Z mu
    sqrtds = zeros(svp,1);
    for i = 1:svp
        if diagS(i) > 1/mu
            sqrtds(i) = sqrt(diagS(i) - 1/mu);
        elseif diagS(i) < -1/mu
            sqrtds(i) = sqrt(diagS(i) + 1/mu);
        else
            sqrtds(i) = 0;
        end
    end
    %sqrtds = sqrt(diagS(1:svp) - 1/mu);
    A.U = A.U(:, 1:svp) * diag(sqrtds);
    A.V = A.V(:, 1:svp) * diag(sqrtds);
    
    X = A.U*A.V';
    E = D - X - 1/mu*Y_1;
    fro_sum = 0;
    for i = 1:p
        fro_sum = fro_sum + E(I(i),J(i))^2;
    end
    if sqrt(fro_sum) > noise_tol
        for i = 1:p
            E(I(i),J(i)) = noise_tol/sqrt(fro_sum) * E(I(i),J(i));
        end
    end  
        
    residual = X + E - D;
    Y_1 = Y_1 + mu*residual;

    %% stop Criterion    
    stopCriterion = norm(residual, 'fro')/d_norm;
    if stopCriterion < tol
        converged = true;
    end    

    %% Maximum iterations reached
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
    
    %% update mu
%     F = sparse(I,J,reshape(matrix_projection(A.U*A.V', observed_entries),[p,1]),m,n) - A.U*A.V';
%     if min(mu, sqrt(mu))*norm(F - E,'fro')/d_norm < 1e-7
%         mu = mu*rho;
%     end
%     E = F;
    mu = mu*rho;
    
    primal_residual(iter) = stopCriterion;
    error(iter) = (norm(X - X_org, 'fro')/norm(X_org, 'fro'))^2;
    muk(iter+1) = mu;
    rak(iter) = rank(X);
end





  