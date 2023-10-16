function demo_cg(S,a0,K,w0,tol)
    if nargin < 5
        tol = 1e-20;
    end
    if nargin < 4
        w0 = zeros(length(a0),1);
    end
    if nargin < 3
        K = 20;
    end
    err_arr = zeros(1,K);
    exact = S\a0;
    k = 1;
    while k <= K
        w_hat = conjugate_gradient_method(w0,S,a0,1e-16,k);
        err_arr(k) = norm(S*w_hat-a0,2)/length(a0);
        k = k + 1;
    end
    figure;
    semilogy(linspace(1,K,K),err_arr);
    grid on
    xlabel('Krylov Subspace Dimension')
    ylabel('||S\hat{w}-a0||_2')
end