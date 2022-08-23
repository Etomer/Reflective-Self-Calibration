function [xopt,yopt,res,jac]=bundletoa2(D,I,J,xt,yt,debug);
% [xopt,yopt]=bundletoa(D,I,J,xt,yt,debug);
%

if nargin<6,
  debug = 0;
end;

for kkk = 1:30;
%     fprintf('Iteration: %d\n',kkk);
    %kkk
    %keyboard;
    [res,jac]=calcresandjac(D,I,J,xt,yt);
    %dz = -(jac\res);
    %dz = -(jac'*jac+speye(size(jac,2)))\(jac'*res);
    dz = -(jac'*jac+0.001*speye(size(jac,2)))\(jac'*res);
    %[u,s,v]=svd(full(jac),0);
    %keyboard;
    %nrparam = size(jac,2);
    %dof = nrparam-6;
    %u = u(:,1:dof);
    %s = s(1:dof,1:dof);
    %v = v(:,1:dof);
    %dz = -v*inv(s)*u'*res;
    [xtn,ytn]=updatexy(xt,yt,dz);
    [res2,jac2]=calcresandjac(D,I,J,xtn,ytn);
    aa = [norm(res) norm(res+jac*dz) norm(res2)];
    bb = aa;
    bb=bb-bb(2);
    bb = bb/bb(1);
    cc = norm(jac*dz)/norm(res);
    % Check that the error actually gets smaller
    if norm(res)<norm(res2),
        % bad
        % check that the update is sufficiently big
        % otherwise it is maybe just numerical limitations
        if cc>1e-4,
            % OK then it is probably just the linearization that
            % is not good enough for this large step size, decrease
            kkkk = 1;
            while (kkkk<50) & (norm(res)<norm(res2)),
                dz = dz/2;
                [xtn,ytn]=updatexy(xt,yt,dz);
                [res2,jac2]=calcresandjac(D,I,J,xtn,ytn);
                kkkk = kkkk+1;
            end
        end
    end
    if debug,
        aa = [norm(res) norm(res+jac*dz) norm(res2)];
        bb = aa;
        bb=bb-bb(2);
        bb = bb/bb(1);
        cc = norm(jac*dz)/norm(res);
        %keyboard;
        aa
        bb
        cc
    end;
    if norm(res2)<norm(res)
        xt = xtn;
        yt = ytn;
    else
        %disp([num2str(kkk) '  stalled']);
    end
    
%     if norm(res)-norm(res2) < 1e-6
%         fprintf('Finished optimization after %d iterations\n',kkk);
%         break;
%     end
end;

xopt = xt;
yopt = yt;

function [res,jac]=calcresandjac(D,I,J,x,y);

nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
idd = 1./dd;
res = dd-D;
II = (1:length(I))';
JJ1 = (I-1)*2+1;
JJ2 = (I-1)*2+2;
%JJ3 = (I-1)*3+3;
JJ4 = (J-1)*2+1+2*m;
JJ5 = (J-1)*2+2+2*m;
%JJ6 = (J-1)*3+3+3*m;

VV1 = idd.*Vt(:,1);
VV2 = idd.*Vt(:,2);
%VV3 = idd.*Vt(:,3);
VV4 = -idd.*Vt(:,1);
VV5 = -idd.*Vt(:,2);
%VV6 = -idd.*Vt(:,3);

jac = sparse([II;II;II;II],[JJ1;JJ2;JJ4;JJ5],[VV1;VV2;VV4;VV5],nn,2*m+2*n);


function [xny,yny]=updatexy(x,y,dz);

m = size(x,2);
n = size(y,2);
dz1 = dz(1:(2*m));
dz2 = dz((2*m+1):end);
xny = x + reshape(dz1,2,m);
yny = y + reshape(dz2,2,n);



