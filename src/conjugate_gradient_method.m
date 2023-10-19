% function x = conjugate_gradient_method(A,b,tol,K)
% % CONJGRAD  Conjugate Gradient Method.
% %   X = CONJGRAD(A,B) attemps to solve the system of linear equations A*X=B
% %   for X. The N-by-N coefficient matrix A must be symmetric and the right
% %   hand side column vector B must have length N.
% %
% %   X = CONJGRAD(A,B,TOL) specifies the tolerance of the method. The
% %   default is 1e-10.
% %
% % Example (highlight lines between %{ and %}, then press F9):
% %{
%   n = 6000;
%   m = 8000;
%   A = randn(n,m);
%   A = A * A';
%   b = randn(n,1);
%   tic, x = conjgrad(A,b); toc
%   norm(A*x-b)
% %}
% % By Yi Cao at Cranfield University, 18 December 2008
% % Updated on 6 Feb 2014.
% %
%     if nargin < 4
%         K = 20;
%     end
%     if nargin<3
%         tol=1e-10;
%     end
%     x = b;
%     r = b - A*x;
%     if norm(r) < tol
%         return
%     end
%     y = -r;
%     z = A*y;
%     s = y'*z;
%     t = (r'*y)/s;
%     x = x + t*y;
% 
%     for k = 1:K;
%        r = r - t*z;
%        if( norm(r) < tol )
%             return;
%        end
%        B = (r'*z)/s;
%        y = -r + B*y;
%        z = A*y;
%        s = y'*z;
%        t = (r'*y)/s;
%        x = x + t*y;
%     end
%  end
% 
% 
% 
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