% Demonstrates the offset-solver and its limitations
%% Initialization
addpath("../src")

%% The solver works on noiseless data
[z_los,z_nlos,offsets] = setup_tdoa_problem();
est_offsets = tdoa_offset_mirror(z_los,z_nlos);

res = est_offsets - offsets


%% Demonstration of the solvers sensitivity to noise

iters = 1000;
noise_stds = logspace(-7,-2);


reses = zeros(length(noise_stds),iters);

for i = 1:length(noise_stds)
    noise_std = noise_stds(i);
    for j = 1:iters
        [z_los,z_nlos,offsets] = setup_tdoa_problem(noise_std);
        est_offsets = tdoa_offset_mirror(z_los,z_nlos);
    
        res = mean(abs(est_offsets - offsets));
        reses(i,j) = res;
    end
end
%
res_mean = median(reses,2);

loglog(noise_stds,res_mean)
xlabel("Measurement noise")
ylabel("Offset estimation error")




%% Helper functions
function [z_los,z_nlos,offsets] = setup_tdoa_problem(noise_std)
    % Solving for 2 senders and 3 receivers since if we had outliers (we do
    % not in this demonstration) we want to solve using the minimal data

    if nargin < 1
        noise_std = 0;
    end
    n_senders = 2;
    n_receivers = 3;
    
    pos_senders = rand(3,n_senders);
    pos_receivers = rand(3,n_receivers);
    
    d_los = zeros(n_receivers,n_senders); %Distance Line Of Sight
    d_nlos = zeros(n_receivers,n_senders); % Distance Non-Line Of Sight
    
    for index_sender = 1:n_senders
        for index_receiver = 1:n_receivers
            d_los(index_receiver,index_sender) = vecnorm(pos_senders(:,index_sender) - pos_receivers(:,index_receiver));
            d_nlos(index_receiver,index_sender) = vecnorm(pos_senders(:,index_sender) - [1;1;-1].*pos_receivers(:,index_receiver));
        end
    end

    %Adding noise
    d_los = d_los + noise_std*randn(n_receivers,n_senders);
    d_nlos = d_nlos + noise_std*randn(n_receivers,n_senders);
    
    offsets = rand(1,n_senders);
    
    z_los = d_los + ones(n_receivers,1).*offsets;
    z_nlos = d_nlos + ones(n_receivers,1).*offsets;
end