function sols = toa_34_mirror_46(d_los,d_nlos)
% TOA_34_MIRROR_46 Solve minimal 3D ToA mirroring problem
%   sols = TOA_34_MIRROR_46(d_los,d_nlos) solves the (3,4) 3D ToA problem
%       where direct distance measurements are given by d_los and the
%       measurements given in d_nlos are assumed to have bounced of the
%       xy-plane. The solutions are returned as an array of structures
%       containing the estimated receiver and sender positions together
%       with the residual norm. This function is a wrapper for the minimal
%       (4,6) 3D ToA problem.

    d = [d_los; d_nlos];
    d2 = d.^2;

    [rrc,ssc] = toa_46(d2,false);
    
    % TODO: Align solutions with xy-plane.
    % TODO: Verify mirror properties of receivers.
    
    % Evaluate residuals and sort solutions.
    sols = struct('r',rrc,'s',ssc,'res',nan(length(rrc),1));
    ress = zeros(length(sols),1);
    for i = 1:length(sols)
        sol = sols(i);
        est_d = pdist2(sol.r',sol.s');
        ress(i) = norm(est_d(:)-d(:));
        sols(i).res = ress(i);
        
        % Keep only non-mirrored receivers.
        sols(i).r = sol.r(:,1:3);
    end
    [~,inds] = sort(ress);
    sols = sols(inds);
end

