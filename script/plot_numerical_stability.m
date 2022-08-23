%% Initialization.
iters = 10000; %Adjust for less noise in graph
addpath("../src")


%% Define solvers.
solvers = [];

% solvers(end+1).name = 'TOA (2,2) mirror';
% solvers(end).solv = @toa_22_mirror;
% solvers(end).m = 2;
% solvers(end).n = 2;
% solvers(end).dim = 2;

solvers(end+1).name = 'TOA (3,4) mirror (192 × 206)';
solvers(end).solv = @(d_los,d_nlos) toa_34_mirror(d_los,d_nlos,@solver_toa_34_mirror_r_dist);
solvers(end).m = 3;
solvers(end).n = 4;
solvers(end).dim = 3;

solvers(end+1).name = 'TOA (3,4) mirror (88 × 102)';
solvers(end).solv = @(d_los,d_nlos) toa_34_mirror(d_los,d_nlos,@solver_toa_34_mirror_r_dist_all_eqs);
solvers(end).m = 3;
solvers(end).n = 4;
solvers(end).dim = 3;

solvers(end+1).name = 'TOA (6,4)';
solvers(end).solv = @toa_34_mirror_46;
solvers(end).m = 3;
solvers(end).n = 4;
solvers(end).dim = 3;


%% Run tests.
% Disable warnings when measuring execution time. This is restored below.
% w = warning('off','all');

ress = inf(length(solvers), iters);
exectime = zeros(1, length(solvers));
for isolv = 1:length(solvers)
    solver = solvers(isolv);
    m = solver.m;
    n = solver.n;
    dim = solver.dim;
    
    RM = eye(dim);
    RM(end) = -1;

    fprintf('(%d/%d) Evalutating %s...\n', isolv, length(solvers), solver.name);

    totaltime = 0;
    for iter = 1:iters
        % Setup problem.
        r = randn(dim, m);
        s = randn(dim, n);

%         r(3,:) = r(3,:)+4;
%         s(3,:) = s(3,:)+4;

        r(rot90(triu(true(m,dim),1))) = 0;
        rr = RM*r;

        d_los = pdist2(r', s');
        d_nlos = pdist2(rr', s');
    
        % Run solver and measure execution time.
        tstart = tic();
        sols = solver.solv(d_los, d_nlos);
        totaltime = totaltime + toc(tstart);

        if isempty(sols)
            continue;
        end

        ress(isolv, iter) = sols(1).res;
    end

    exectime(isolv) = totaltime / iters;
end

% Restore warnings.
% warning(w);

%% Plot numerical stability.
figure(1);
edges = -15:0.5:3;
% edges = -17:0.5:3;
centers = edges(1:end-1) + diff(edges) / 2;
binwidth = edges(2)-edges(1);
for isolv = 1:length(solvers)
    solver = solvers(isolv);

    hc = histcounts(log10(ress(isolv, :)), edges) / iters / binwidth;
    plot(centers, hc, 'LineWidth', 2);
    hold on

    fprintf('%s failed to find a real solution in %.1f %% of all tests.\n', ...
        solver.name, 100*sum(~isfinite(ress(isolv, :)))/iters);
end
hold off

legend(solvers.name, 'Interpreter', 'none');
xlabel('log_{10}(error)');
xlim(edges([1 end]));
ylim([0 0.2]);
ylabel('Relative Frequency');

set(gca, 'FontName', 'Times');
set(gca, 'FontSize', 12);

%% Plot execution time.
figure(2);
bp = bar(1000*exectime, 'FaceColor', 'flat');
bp.CData = lines(length(solvers));
ylabel('Execution time [ms]');
xticklabels({solvers.name});
xtickangle(45);
set(gca, 'TickLabelInterpreter', 'none');
title('Execution time');

fprintf('Execution times (microseconds):\n');
fprintf('%5.0f\n', 1e6*exectime);
