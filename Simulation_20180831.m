close all;

tic
%%%%%%%%%%%%%%%%
% Parameters   %
%%%%%%%%%%%%%%%%

num_Step = 10;
iter = 1000;
Type =2;
res = zeros(iter,num_Step);
maxX = 100;
Wifi_range = 20;


%% simulation start
for ii = 1:iter
    for N = num_Step:10:10*num_Step
        Nloc =  rand(N, 2) * maxX;
        for n = 1:N
            if(Nloc(n,1)<Wifi_range && Nloc(n,2)<Wifi_range)
                res(ii, N/10) = res(ii, N/10) + 1;
            else
                res(ii, N/10) = res(ii, N/10) + 1/N;
            end
            
        end
        res(ii,N/10) = res(ii,N/10) / N;
    end
end

final_res = mean(res);


toc

plot(final_res)


% final_res = total_res / iter;
% throughput = final_res*Trans_slot*aSlotTime/total_time;
% 
% x = 100:100:num_GR*100;
% x = 1./x *data_time/aSlotTime; 
% figure(1);
% plot(x, throughput(1,:),'-.bo');
% hold on;
% plot(x, throughput(2,:), '-.r+');
% hold on;
% plot(x, throughput(3,:), '-.g*');
% xlabel('Packet Generation rate'); ylabel('Normalized System Throughput');
% legend('Uniform distribution based conventional scheme','Linear distribution based proposed scheme', 'Exponential distribution based proposed scheme');

% x = 1:num_STA;
% figure(1);
% plot(x, success(:),'-.b*');
% xlabel('Node Index'); ylabel('Number of Successful transmission');