function sols = toa_34_mirror_fromheights(dp,z,w)
% TOA_34_MIRROR Solve minimal 3D ToA mirroring problem
%   sols = TOA_34_MIRROR(d_los,d_nlos) solves the (3,4) 3D ToA problem
%       where direct distance measurements are given by d_los and the
%       measurements given in d_nlos are assumed to have bounced of the
%       xy-plane. The solutions are returned as an array of structures
%       containing the estimated receiver and sender positions together
%       with the residual norm.

%     dm = (d_nlos.^2-d_los.^2)/4;
%     dp = (d_nlos.^2+d_los.^2)/2;
% 
%     [U,S,V] = svd(dm);
%     z = -sqrt(S(1))*U(:,1);
%     w = -sqrt(S(1))*V(:,1);

    data = [dp(:); z.^2; w.^2];
    sols_alpha = solver_toa_34_mirror_r_dist(data);
%     sols_alpha = solver_toa_34_mirror_r_dist_all_eqs(data);
    keep = max(abs(imag(sols_alpha)./real(sols_alpha))) < 1e-6;
    keep = keep & all(real(sols_alpha) >= 0);
    good_sols = real(sols_alpha(:,keep));

    sols = struct('r',{},'s',{},'res',{},'beta',{});
    for i = 1:size(good_sols,2)
        dr2 = good_sols(1:3,i);
        alpha = good_sols(4,i);
        beta = sqrt(alpha);
        
        r = zeros(3,3);
        s = zeros(3,4);
        
        % Find heights.
        r(3,:) = z*beta;
        s(3,:) = w/beta;
        
        % Find receivers in 2D.
        x1 = sqrt(dr2(3));
        x3 = (dr2(3)-dr2(1)+dr2(2))./(2*x1);
        x2 = sqrt(dr2(2)-x3.^2);
        if ~isreal(x1) || ~isreal(x2)
            continue;
        end
        r(2,2) = x1;
        r(1,3) = x2;
        r(2,3) = x3;
        
        % Trilaterate sender positions in 2D.
        d_2d = dp-r(3,:)'.^2-s(3,:).^2;
        if any(d_2d < 0,'all')
            continue;
        end
        Rt = r(1:2,2:3);
        Dt = d_2d(2:3,:)-d_2d(1,:);
        A = 2*Rt';
        b = dot(Rt,Rt)'-Dt;
        s(1:2,:) = A\b;
        
        sols(end+1) = struct('r',r,'s',s,'res',nan,'beta',beta);
    end
    
    % Evaluate residuals and sort solutions.
%     ress = zeros(length(sols),1);
%     for i = 1:length(sols)
%         sol = sols(i);
%         ss = diag([1 1 -1])*sol.s;
%         est_d = pdist2(sol.r',[sol.s ss]');
%         ress(i) = norm(est_d(:)-[d_los(:); d_nlos(:)]);
%         sols(i).res = ress(i);
%     end
%     [~,inds] = sort(ress);
%     sols = sols(inds);
end

