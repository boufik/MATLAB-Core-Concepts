function [beta_dual, b_dual, l] = svm_dual(X, y, P, q, lambda)
    
    [n, p] = size(X);
    cvx_begin quiet
        variable l(n);
        minimize(0.5 * l'*P*l - q'*l);
        subject to
        l >= 0;
        y'*l == 0;
    cvx_end
    
    beta_dual = zeros(1, p);
    b_dual = 0;
    for k = 1 : n
        beta_dual = beta_dual + lambda(k) * X(k, :) * y(k);
    end
    for k = 1 : n
        b_dual = b_dual + (1/n) * (y(k) - beta_dual*X(k, :)');
    end

end

