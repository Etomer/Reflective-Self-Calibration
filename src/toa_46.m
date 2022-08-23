function [rrc,ssc] = toa_46(d2,refine)
% TOA46 Solves the minimal (4,6) time of arrival problem.
%   [rrc,ssc] = TOA_46(d2,refine) solves the minimal (4,6) ToA problem
%       given a matrix d2 of squared distances. d2 can be of size [4 6] or
%       [6 4] depending on the number of receivers and six senders. rrc and
%       ssc are cell arrays containing the solutions for receiver and
%       sender positions. refine specifies whether or not the solutions
%       should be refined by local optimization.
%
%   For details please see:
%       V. Larsson (2017) Polynomial Solvers for Saturated Ideals
%       Y. Kuang (2013) A complete characterization and solution to the
%           microphone position self-calibration problem

    if isequal(size(d2),[4 6])
        % Solver expects 6 receivers and 4 senders.
        d2 = d2';
        swapRS = true;
    elseif isequal(size(d2),[6 4])
        swapRS = false;
    else
        error('Expected matrix of size [4 6] or [6 4].');
    end
    
    if nargin < 2
        refine = false;
    end

    cl = compactionmatrix(6);
    cr = compactionmatrix(4);

    B = -cl*d2*cr'/2;

    % Factor B in R and S.
    R = B';
    S = eye(3);

    da = d2(1,1);
    db = cr*d2(1,:)';
    dc = cl*d2(:,1);

    % Find linear requirements on H and b.
    indsI = [1 2 3 1 1 2];
    indsJ = [1 2 3 2 3 3];
    mult = [1 1 1 2 2 2]; % Due to symmetry in H.
    
    AA = [mult.*R(indsI,:)'.*R(indsJ,:)' -2*R'];
    bb = dc;

    Hbc = AA\bb;
    [~,~,v] = svd(AA);
    Hbb = v(:,6:9);

    % Solve polynomial system using action matrix method.
    data = [Hbb(:); Hbc; db; da];
    sols = solver_toa_46(data);

    rrc = cell(1,size(sols,2));
    ssc = cell(1,size(sols,2));
    for i=1:size(sols,2)
        x = sols(1:4,i);
        
        if ~isreal(x)
            continue;
        end
        
        if refine
            x = refine_toa_46(Hbb,Hbc,db,da,x,10,1e-6);
        end
        
        Hb = Hbb*x+Hbc;
        
        H = Hb([1 4 5; 4 2 6; 5 6 3]);
        b = Hb([7; 8; 9]);

        try
            Lti = chol(H);
        catch
            continue;
        end
        
        rrc{i} = Lti*[zeros(3,1) R];
        ssc{i} = Lti'\([zeros(3,1) S]+repmat(b,1,4));
    end
    rrc(cellfun(@isempty,rrc)) = [];
    ssc(cellfun(@isempty,ssc)) = [];
    
    if swapRS
        temp = rrc;
        rrc = ssc;
        ssc = temp;
    end
end

function [cc,dd] = compactionmatrix(n)
    cc = [-ones(n-1,1) eye(n-1)];
    dd = [1 zeros(1,n-1);cc];
end

