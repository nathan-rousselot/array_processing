function w = conjugate_gradient_method(w0, S, a0, tol, max_iter)
    if nargin < 5
        max_iter = 1000;
    end
    if nargin < 4
        tol = 1e-10;
    end
    [alpha,beta] = deal(zeros(1,max_iter));
    [r,w,p] = deal(ones(length(w0),max_iter));
    w(:,1) = w0;
    r(:,1) = a0-S*w(:,1);
    p(:,1) = r(:,1);
    if norm(r(:,1),2) <= tol
        w = w(:,1);
        return;
    end
    k = 1;
    while (k <= max_iter)
        alpha(k) = (r(:,k)'*r(:,k))/(p(:,k)'*S*p(:,k));
        w(:,k+1) = w(:,k)+alpha(k)*p(:,k);
        r(:,k+1) = r(:,k)-alpha(k)*S*p(:,k);
        if norm(r(:,k+1),2) <= tol
            w = w(:,k+1);
            return;
        end
        beta(k) = (r(:,k+1)'*r(:,k+1))/(r(:,k)'*r(:,k));
        p(:,k+1) = r(:,k+1)+beta(k)*p(:,k);
        k = k + 1;
    end
    w = w(:,k);
end