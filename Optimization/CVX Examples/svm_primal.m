function [beta_primal, b_primal] = svm_primal(X, y)

    [n, p] = size(X);
    cvx_begin quiet
        variables beta_primal(p) b_primal;
        minimize(0.5 * norm(beta_primal))
        subject to
            y .* (X*beta_primal + b_primal) >= 1;
    cvx_end

end

