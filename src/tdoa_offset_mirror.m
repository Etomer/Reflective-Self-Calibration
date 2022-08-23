function o = tdoa_offset_mirror(z_los,z_nlos,refine,maxiter)
% TDOA_OFFSET_MIRROR Find TDoA mirror offsets
%   o = TDOA_OFFSET_MIRROR(z_los,z_nlos) finds the TDoA offsets in the
%       measurements such that the distance measurements are given by
%       d_los = z_los-o, d_nlos = z_nlos-o. The solver utilizes the fact
%       that dm = d_nlos.^2-d_los.^2 has rank 1.
%   o = TDOA_OFFSET_MIRROR(z_los,z_nlos,refine,maxiter) enables iterative
%       refinment of the solution using maxiter (default 10) iterations of
%       coordinate descent.

    [m,n] = size(z_los);
    
    if nargin < 3
        refine = false;
    end
    if nargin < 4
        maxiter = 10;
    end
    
    % dm = 0.25*((z_nlos-o).^2-(z_los-o).^2) has rank 1 and is linear in o.
    dmc = 0.25*(z_nlos.^2-z_los.^2); % Constant term.
    dmo = -0.5*(z_nlos-z_los); % Offset coefficient.
    
    permi = nchoosek(1:m,2);
    permj = nchoosek(1:n,2);
    
    nmon = n*(n-1)/2+n;
    nmon2 = n*(n-1)/2;
    neqs = size(permi,1)*size(permj,1);
    assert(neqs >= nmon,...
        'Linear system is underdetermined for %d receivers and %d senders',m,n);
    A = zeros(neqs,nmon);
    b = zeros(neqs,1);
    
    % dm has rank 1 <=> all minors of order 2 are zero
    % Each minor has four monomials: o_j1*o_j2, o_j1, o_j2, and 1, for some
    % indices [j1 j2] = cs.
    k = 1;
    for i = 1:size(permi,1)
        rs = permi(i,:);
        for j = 1:size(permj,1)
            cs = permj(j,:);
            A(k,j) = det(dmo(rs,cs));
            A(k,nmon2+cs(1)) = det([dmo(rs,cs(1)) dmc(rs,cs(2))]);
            A(k,nmon2+cs(2)) = det([dmc(rs,cs(1)) dmo(rs,cs(2))]);
            b(k) = -det(dmc(rs,cs));
            k = k+1;
        end
    end
    
    % A*mv = b where mv contains second and first order monomials in o.
    mv = A\b;
    o = mv(nmon2+(1:n))';

    % Refine solution using coordinate descent.
    % We seek to minimizer norm(dmc+dmo*diag(o)-g*h','fro').
    if refine
        % Find initial g and h.
        dm = dmc+dmo*diag(o);
        [U,S,V] = svd(dm);
        g = sqrt(S(1,1))*U(:,1);
        h = sqrt(S(1,1))*V(:,1);
        
        for iter = 1:maxiter
            % Solve for o and h.
            for j = 1:n
               oh = -([dmo(:,j) -g] \ dmc(:,j));
               o(j) = oh(1);
               h(j) = oh(2);
            end
            
            % Solve for o and g.
            A = zeros(m*n,n+m);
            for j = 1:n
                A(m*(j-1)+1:m*j,j) = dmo(:,j);
                A(m*(j-1)+1:m*j,n+1:end) = -h(j)*eye(m);
            end
            og = -(A \ dmc(:));
            o = og(1:n)';
            g = og(n+1:end);
            
            % Adjust the norm of g and h.
            adjust = sqrt(norm(g)/norm(h));
            g = g ./ adjust;
            h = h .* adjust;
        end
    end
end
