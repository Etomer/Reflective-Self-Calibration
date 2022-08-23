function [R,sols,res] = solver_opttrilat(s0,d,w,use_eigenvector)
% Finds all stationary points to
% min_x sum_j w(j) * (|x-s(:,k)|^2 - d(k)^2)^2

if nargin < 4
    use_eigenvector = 1;
end

n = size(s0,1);
d = d(:);
if nargin < 3 || isempty(w)
    w = 1./d.^2;
else
    w = w(:);
end
% Normalize weights
w = w / sum(w);


% shift center
t = sum(s0 .* w(:,ones(1,n))',2);
s = s0 - t * ones(1,size(s0,2));

% Construct A,b such that (x'*x)*x + A*x + b = 0
ws2md2 = w.*(sum(s.^2,1).'-d.^2);
A = 2*(s .* w(:,ones(n,1)).') * s.' + sum(ws2md2) * eye(n);
b = -s*ws2md2;

[V,D] = eig(A);
bb = V'*b;

% basis = [x^2,y^2,z^2,x,y,z,1]
AM = [-D diag(-bb) zeros(n,1);
     zeros(n,n) -D -bb;
     ones(1,n) zeros(1,n+1)];

if use_eigenvector
    [VV,DD] = eig(AM);
    VV = VV ./ (ones(2*n+1,1)*VV(end,:));
    r = V*VV(n+1:2*n,:);
else
    DD = eig(AM);
    % eigenvector-less solution extraction
    r = zeros(n,2*n+1);
    for k = 1:2*n+1
        z = [zeros(n,n);-eye(n)];
        T = AM - DD(k)*eye(2*n+1);
        r(:,k) = (T(:,1:end-1).' \ z).'*T(:,end);
    end
    r = V*r;
end


% perform some refinement on the roots
for i = 1:2*n+1    
    ri = r(:,i);
    if max(abs(imag(ri))) > 1e-6
        continue
    end
    for k = 1:3
        res = (ri'*ri)*ri + A*ri + b;
        if norm(res) < 1e-8
            break;
        end
        J = (ri'*ri)*eye(n) + 2 *(ri*ri') + A;
        ri = ri - J\res;
    end
    r(:,i) = ri;
end

% Revert translation of coordinate system
sols = r + t*ones(1,size(r,2));

% find best stationary point
cost = inf;
R = zeros(n,1);
for k = 1:size(sols,2)
    if sum(abs(imag(sols(:,k)))) > 1e-6
        continue;
    end
    rk = real(sols(:,k));
    cost_k = sum(w'.*(sum((rk(:,ones(1,size(s0,2)))-s0).^2)- d'.^2).^2);
    if cost_k < cost
        cost = cost_k;
        R = rk;
    end
end

if nargout == 3
    lambda = diag(DD);
    res = [];
    for k = 1:2*n+1
        rk = r(:,k);
        lambdak = lambda(k);
        res(:,k) = [(rk.'*rk)*rk + A*rk + b; lambdak - rk.'*rk];
    end
end


