%% load and setup data

addpath("../src/")
addpath("../src/solvers/")
addpath("../data/")
load my_0014_data_v2

M = size(rgt,2);
N = size(sgt_resamp,2);

 bnd = 0.05;
%bnd = 0.1; % use tracks that are within bnd of groundtruth

anchor = 6;
uu = u_meters(anchor,:); % use 6th mic as offset

rgtm = rgt;
rgtm(end,:)=-rgt(end,:); % gt mirror positions

dd6 = sqrt(sum((rgt(:,anchor)-sgt_resamp).^2)); % distance from sender to mic 6

uu2 = uu; % move all measurements with offset
for iii = 1:M
    uu2{iii} = uu{iii}+dd6;
end


ddgt = uu; % ground truth distances from sender to mics
ddmgt = uu; % ground truth distances from sender to mirror mics

for iii = 1:M
    ddgt{iii} = sqrt(sum((rgt(:,iii)-sgt_resamp).^2));
    ddmgt{iii} = sqrt(sum((rgtm(:,iii)-sgt_resamp).^2));    
end

uu_clean = uu2; % remove outliers based on gt
for iii = 1:12
    d1 = ddgt{iii};
    d2 = ddmgt{iii};
    for jjj = 1:4
        tmp = uu2{iii}(jjj,:);
        ok = abs(tmp-d1)<bnd | abs(tmp-d2)<bnd;
        tmp(~ok) = nan;
        uu_clean{iii}(jjj,:)=tmp;
    end
end



D = nan(M,N);
dd_sub = nan(M,N);
ddm_sub = nan(M,N);
for iii = 1:M
    tmp = sort(uu_clean{iii});
    ok1 = sum(isfinite(tmp))==2;
    ok2 = abs(tmp(2,:)-tmp(1,:))>2*bnd;
    ok3 = abs(tmp(3,:)-tmp(1,:))>2*bnd;
    ok4 = sum(isfinite(tmp))==3;
    D(iii,ok1 & ok2) = (tmp(2,ok1 & ok2).^2-tmp(1,ok1 & ok2).^2)/4;
    D(iii,ok3 & ok4) = (tmp(3,ok3 & ok4).^2-tmp(1,ok3 & ok4).^2)/4;
    dd_sub(iii,ok1 & ok2)=tmp(1,ok1 & ok2);
    ddm_sub(iii,ok1 & ok2)=tmp(2,ok1 & ok2);
    dd_sub(iii,ok3 & ok4)=tmp(1,ok3 & ok4);
    ddm_sub(iii,ok3 & ok4)=tmp(3,ok3 & ok4);   
end



Dgt = zeros(M,N);
for iii = 1:M
    Dgt(iii,:)=(ddmgt{iii}.^2-ddgt{iii}.^2)/4;
end


X = D(setdiff(1:M,anchor),:);
Xgt = Dgt(setdiff(1:M,anchor),:);
  
ok = sum(isfinite(X))>2;
X = X(:,ok);
Xgt = Xgt(:,ok);
W = isfinite(X);

[M2,N2] = size(X);

dd_sub = dd_sub(setdiff(1:M,anchor),ok);
ddm_sub = ddm_sub(setdiff(1:M,anchor),ok);

rgt2 = rgt;
rgt2(:,anchor) = [];
sgt2 = sgt_resamp(:,ok);


%% run low rank estimation, uses https://github.com/hamburgerlady/miss-ranko


params.rk = 1;
params.robust = 0;
params.inlierbnd= 0.5000;
params.nrinliersbnd= 1;
params.initnn= 1;
params.minhalfn= 1;
params.bundleiter= 2;
params.glueraniter= 5;
params.extendraniter= 5;
params.maxiter= 5000;
params.finalnormbnd= 100;
params.cutty= 0.7500;
params.gksigge= 3;
params.gksize= 20;
params.maxstatic= 1;
params.finN= M2;
params.finM= N2;

sol1 = mr_solver_rankn(X,W,params);
sol2 = mr_solver_rankn(Xgt,W,params);

Ue = sol1.U(:);
Ve = sol1.V;
Ve = Ve(:);

Ugt = rgt(3,setdiff(1:M,anchor));
Ugt = Ugt(:);
Vgt = sgt_resamp(3,ok);
Vgt = Vgt(:);

s1 = mean(Vgt./Ve); 
s2 = mean(Ue./Ugt); 
sc = mean([s1 s2]);


figure(1);
clf
ll = plot(Ve*sc);
hold on
set(ll,'LineWidth',2);
ll = plot(Vgt,'--');
set(ll,'LineWidth',2);
legend({'Estimated','Ground Truth'});
ylabel('Height (meters)');
xlabel('Sender');
axis([1 sum(ok) 0 2])

figure(2)
clf
plot(Ue/sc,'*')
hold on
plot(Ugt,'o')
legend({'Estimated','Ground Truth'});
ylabel('Height (meters)');
xlabel('Receiver');
axis([1 M 0 2])


%% initialize a large number of solutions for 2D estimation and low-rank scale
Ds = nan(M,N);

for iii = 1:M
    tmp = sort(uu_clean{iii});
    ok1 = sum(isfinite(tmp))==2;
    ok2 = abs(tmp(2,:)-tmp(1,:))>2*bnd;
    ok3 = abs(tmp(3,:)-tmp(1,:))>2*bnd;
    ok4 = sum(isfinite(tmp))==3;
    Ds(iii,ok1 & ok2) = (tmp(2,ok1 & ok2).^2+tmp(1,ok1 & ok2).^2)/2;
    Ds(iii,ok3 & ok4) = (tmp(3,ok3 & ok4).^2+tmp(1,ok3 & ok4).^2)/2;
end

Dsgt = zeros(M,N);
for iii = 1:M
    Dsgt(iii,:)=(ddmgt{iii}.^2+ddgt{iii}.^2)/2;
end

Ds = Ds(setdiff(1:M,anchor),ok);
Dsgt = Dsgt(setdiff(1:M,anchor),ok);

allsols = struct('r',{},'s',{},'res',{},'beta',{},'idi',{},'idj',{});

niter = 2000;


for iter = 1:niter
    

    idj0 = randi(349);
    idi = find(isfinite(Ds(:,idj0)));
    idi = idi(randperm(length(idi),3));
    idj = find(sum(isfinite(Ds(idi,:)))==3);
    if length(idj)>3
        idj = idj(randperm(length(idj),4));
        
      
        Di = Ds(idi,idj);
        
        if min(svd(Di))>1
        zi = Ue(idi);
        wi = Ve(idj);
    
        sols = toa_34_mirror_fromheights(Di,zi,wi);
        for si = 1:length(sols)
            sols(si).idi = idi;
            sols(si).idj = idj;
            allsols = [allsols;sols(si)];
        end
        end
        
    end
    
    
end

nsol = length(allsols);
allbetas = zeros(1,nsol);
for iii = 1:nsol
    allbetas(iii) = allsols(iii).beta;
end




%% use voting scheme to finds good scale (beta)
nbins = 10;
betabnd = 0.05;
nsol = length(allsols);
betas = zeros(1,nsol);
maxy = zeros(1,nsol);
for iii = 1:nsol
    betas(iii) = allsols(iii).beta;
    maxy(iii) = max([Ue*allsols(iii).beta;Ve/allsols(iii).beta]);
end

badids = maxy>2;


bb = linspace(min(betas(~badids)),max(betas(~badids)),nbins);
hh = hist(betas(~badids),bb);
[M,ids] = max(hh);
betaopt = bb(ids);

goodids = abs(betas-betaopt)<betabnd;
allsols2 = allsols(goodids);


dd2opt = Ds-(Ue*betaopt).^2-(Ve'/betaopt).^2;

%% expand tentative solutions and choose best one
bestres = inf;
for solind = 1:length(allsols2)
sol = allsols2(solind);

for iter = 1:10
sol.idj = [];
sol.s = [];
for jjj = 1:N2
    dsub = sqrt(dd2opt(sol.idi,jjj));
    oksub = isfinite(dsub);
    if sum(oksub)>2
        sol.idj = [sol.idj jjj];
        se = solver_opttrilat(sol.r(1:2,oksub),dsub(oksub),1./dsub(oksub).^2);
        %disp(sols)
        %disp('hepp')
        %[se,e] = choose_opttrilat(sol.r(1:2,oksub),dsub(oksub));
        sol.s = [sol.s se];
    end
end

sol.idi = [];
sol.r = [];
for iii = 1:M2
    dsub = sqrt(dd2opt(iii,sol.idj));
    oksub = isfinite(dsub);
    if sum(oksub)>2
        sol.idi = [sol.idi iii];
        re = solver_opttrilat(sol.s(1:2,oksub),dsub(oksub),1./dsub(oksub).^2);
        %disp(sols)
        %disp('hepp')
        %[se,e] = choose_opttrilat(sol.r(1:2,oksub),dsub(oksub));
        sol.r = [sol.r re];
    end
end

end

if length(sol.idi)==11 && length(sol.idj)==349
inliers = isfinite(dd2opt);
[I,J,~]=find(inliers);
ind = sub2ind(size(dd2opt),I,J);
D = sqrt(dd2opt(ind));
[xopt,yopt,res,jac]=bundletoa2(D,I,J,sol.r,sol.s);

for iii = 1:4
   [xopt,yopt,res,jac]=bundletoa2(D,I,J,xopt,yopt);
    %disp(res'*res);
end

resi = res'*res;
if resi < bestres
    bestres = resi;
    bestsol = sol;
    bestsol.xopt = xopt;
    bestsol.yopt = yopt;
end


end

end



%% normalise 2D estimate and ground truth
[x,y] = toa_normalise(bestsol.xopt,bestsol.yopt);
[rtmp,stmp] = toa_normalise(rgt2(1:2,:),sgt2(1:2,:));




%% bundle from gt

inliers = isfinite(dd2opt);
[I,J,~]=find(inliers);
ind = sub2ind(size(dd2opt),I,J);
D = sqrt(dd2opt(ind));
[xoptgt,yoptgt,resgt,jac]=bundletoa2(D,I,J,rgt2(1:2,:),sgt2(1:2,:));


%% plot results

fs = 21;
h = figure(1);
clf
hold on
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', fs);

for iii = 1:M2
    ll = plot(dd_sub(iii,:),'o');
    llm = plot(ddm_sub(iii,:),'x');
    set(llm,'color',get(ll,'color'));
end
legend({'Direct distance','Mirror distance'})
xlabel('Sender number');
ylabel('distance (m)')
axis([1 N2 0 5.5])
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'../figs/real_input_data','-dpdf','-r0')



h = figure(2);
clf
hold on
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', fs);

llr = plot(x(1,:),x(2,:),'*');
set(llr,'MarkerSize',10);
set(llr,'LineWidth',1.5);
lls = plot(y(1,:),y(2,:),'x');
set(lls,'MarkerSize',10);
set(lls,'LineWidth',1.5);
set(lls,'color',get(llr,'color'));
llr = plot(rtmp(1,:),rtmp(2,:),'o');
set(llr,'MarkerSize',10);
set(llr,'LineWidth',1.5);
llsgt = plot(stmp(1,:),stmp(2,:),'+');
set(llsgt,'MarkerSize',10);
set(llsgt,'color',get(llr,'color'));
set(llsgt,'LineWidth',1.5);
legend({'Est. receivers','Est. senders','G.T. receivers','G.T. senders'},'Location','NorthWest')
axis equal
xlabel('(m)');
ylabel('(m)');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'result_real_2d','-dpdf','-r0')


h = figure(3);
clf
hold on
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', fs);

llr = plot(Ue/sc,'*');
set(llr,'MarkerSize',10);
set(llr,'LineWidth',1.5);

llgt = plot(Ugt,'o');
set(llgt,'MarkerSize',10);
set(llgt,'LineWidth',1.5);
set(llgt,'color',get(llsgt,'color'));

legend({'Estimated','Ground Truth'});
axis([1 M2 0 2])
ylabel('Height (meters)');
xlabel('Receiver');

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'result_real_height_receiver','-dpdf','-r0')


h = figure(4);
clf
hold on
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', fs);

lls = plot(Ve*sc,'x');
set(lls,'MarkerSize',10);
set(lls,'LineWidth',1.5);
%set(lls,'color',get(llr,'color'));
llgt = plot(Vgt,'+');
set(llgt,'MarkerSize',10);
set(llgt,'color',get(llsgt,'color'));
set(llgt,'LineWidth',1.5);
axis([1 N2 0 2])
legend({'Estimated','Ground Truth'});
ylabel('Height (meters)');
xlabel('Sender');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'result_real_height_sender','-dpdf','-r0')



h = figure(5);
clf
hold on
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', fs);

llr = plot(Ue*betaopt,'*');
set(llr,'MarkerSize',10);
set(llr,'LineWidth',1.5);

llgt = plot(Ugt,'o');
set(llgt,'MarkerSize',10);
set(llgt,'LineWidth',1.5);
set(llgt,'color',get(llsgt,'color'));

legend({'Estimated','Ground Truth'});
axis([1 M2 0 2])
ylabel('Height (meters)');
xlabel('Receiver');

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'result_real_height_receiver_bestbeta','-dpdf','-r0')

h = figure(6);
clf
hold on
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', fs);

lls = plot(Ve/betaopt,'x');
set(lls,'MarkerSize',10);
set(lls,'LineWidth',1.5);
%set(lls,'color',get(llr,'color'));
llgt = plot(Vgt,'+');
set(llgt,'MarkerSize',10);
set(llgt,'color',get(llsgt,'color'));
set(llgt,'LineWidth',1.5);
axis([1 N2 0 2])
legend({'Estimated','Ground Truth'});
ylabel('Height (meters)');
xlabel('Sender');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'result_real_height_sender_bestbeta','-dpdf','-r0')



%% plot 3D estimate and gt
Rgt = [rtmp;Ugt'];
Sgt = [stmp;Vgt'];
mmh = mean([Rgt(3,:) Sgt(3,:)]);
Rgt = Rgt-[0;0;mmh];
Sgt = Sgt-[0;0;mmh];


Re = [x;Ue'*betaopt];
Se = [y;Ve'/betaopt];
mmh = mean([Re(3,:) Se(3,:)]);
Re = Re-[0;0;mmh];
Se = Se-[0;0;mmh];



h = figure(7);
clf
hold on
set(gca, 'FontName', 'Times');
set(gca, 'FontSize', fs);

llr = plot3(Re(1,:),Re(2,:),Re(3,:),'*');
set(llr,'MarkerSize',10);
set(llr,'LineWidth',1.5);
lls = plot3(Se(1,:),Se(2,:),Se(3,:),'x');

set(lls,'MarkerSize',10);
set(lls,'LineWidth',1.5);
set(lls,'color',get(llr,'color'));

llr = plot3(Rgt(1,:),Rgt(2,:),Rgt(3,:),'o');
set(llr,'MarkerSize',10);
set(llr,'LineWidth',1.5);
llsgt = plot3(Sgt(1,:),Sgt(2,:),Sgt(3,:),'+');
set(llsgt,'MarkerSize',10);
set(llsgt,'color',get(llr,'color'));
set(llsgt,'LineWidth',1.5);
legend({'Est. receivers','Est. senders','G.T. receivers','G.T. senders'},'Location','NorthWest')

xlabel('(m)');
ylabel('(m)');
zlabel('(m)');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
view(3)
axis equal

