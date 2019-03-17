close all;

tic
%%%%%%%%%%%%%%%%
% Parameters   %
%%%%%%%%%%%%%%%%

iter = 1000;
Types =2;


maxX = 100;
TransRan = 10;

numAP = 10;
% Nloc = zeros(numN, 3);
% Nloc(:,1:2) = 2*(rand(numN, 2)-0.5) * maxX;
APloc = zeros(numAP, 3);
APloc(:,1:2) = 2*(rand(numAP, 2)-0.5) * maxX;

%% AP relocate
for ap = 1:numAP
    flagDup = 0;
    while(true)
        for apRef = 1:numAP
            if(apRef==ap)
                continue;
            end
            if(sqrt((APloc(apRef,1)-APloc(ap,1))^2 + (APloc(apRef,2)-APloc(ap,2))^2)<TransRan)
                flagDup = 1;
            end
        end
        if(flagDup==0)
            break;
        end
        APloc(ap,1:2) = 2*(rand(1, 2)-0.5) * maxX;
        flagDup = 0;
    end
end

numNStart = 100;
numNStep = 100;
numNMax = 1000;
numNArray = numNStart:numNStep:numNMax;
NN = length(numNArray);

resHard = zeros(iter,NN);
resSoft = zeros(iter,NN);
resNoHet = zeros(iter,NN);


%% simulation start
for ii = 1:iter
    for nn = 1:NN
        numN = numNArray(nn);
        Nloc = zeros(numN, 3);
        Nloc(:,1:2) = 2*(rand(numN, 2)-0.5) * maxX;
        for n = 1:numN
            for ap = 1:numAP
                if (sqrt((Nloc(n,1)-APloc(ap,1))^2 + (Nloc(n,2)-APloc(ap,2))^2)<TransRan)
                    Nloc(n,3) = ap;
                end
            end
        end
        numNResHard = numN - length(find(Nloc(:,3)>0));
        numNResSoft = numN - length(find(Nloc(:,3)>0))/2;
        for n = 1:numN
            if(Nloc(n,3)>0)
                resHard(ii, nn) = resHard(ii, nn) + 1/length(find(Nloc(:,3)==Nloc(n,3)));
                resSoft(ii, nn) = resSoft(ii, nn) + 2/length(find(Nloc(:,3)==Nloc(n,3)));
            else
                resHard(ii, nn) = resHard(ii, nn) + 1/numNResHard;
                resSoft(ii, nn) = resSoft(ii, nn) + 1/numNResSoft;
            end
            resNoHet(ii, nn) = resNoHet(ii, nn) + 1/numN;
        end
        resHard(ii,nn) = resHard(ii,nn) / numN;
        resSoft(ii,nn) = resSoft(ii,nn) / numN;
        resNoHet(ii,nn) = resNoHet(ii,nn) / numN;
    end
end


toc

%% topology plot
figure(1);
plot(0,0,'o');
hold on;
plot(Nloc(:,1), Nloc(:,2),'.');
hold on;
plot(APloc(:,1), APloc(:,2),'x');

figure(2)
final_resHard = mean(resHard);
final_resSoft = mean(resSoft);
final_resNoHet = mean(resNoHet);
plot(final_resHard)
hold on;
plot(final_resSoft)
hold on;
plot(final_resNoHet)

%% results
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