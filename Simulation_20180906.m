close all;
clear variables;

tic
%%%%%%%%%%%%%%%%
% Parameters   %
%%%%%%%%%%%%%%%%

iter = 1;
Types =2;

maxX = 100;
TransRan = 10;

numNStart = 10000;
numNStep = 100;
numNMax = 10000;
numNArray = numNStart:numNStep:numNMax;
NN = length(numNArray);

delayFail = 10;
delaySuc = 200;

resHard = zeros(iter,NN);
resSoft = zeros(iter,NN);
resNoHet = zeros(iter,NN);

%simulation time
maxT = 1000*60*100;
pktGenRate = 1000*60;

% for later use
nn = 1;

%% simulation start
numN = numNArray(nn);
% Register Construction, 1:average Delay, 2: # delivered pkts, 3:#
% packets in a queue, 4:time of first packet, 5: next packet timer, 6: next trial timer
% Register Success: 1st #trials, 4:2nd...
regN = zeros(numN,6);
regSuc = (-1)*ones(numN, 10* maxT/pktGenRate);
GR = ceil(random('exp', ones(numN,1)*pktGenRate));

t = 0;
flagSuc =0;
while (t<maxT)
    
    t= t+1;
    % packet generation
    regN(:,5) = regN(:,5)+ (regN(:,5)==0).* GR;
    regN(:,3) = regN(:,3) + (regN(:,5)==1);
    regN(:,4) = regN(:,4) + t*((regN(:,5)==1)&(regN(:,3)==1));
    regN(:,5) = (regN(:,5)>0).* (regN(:,5)-1);
    
    regN(:,6) = (regN(:,6)>0).* (regN(:,6)-1);
    
    % success timer setup
    if flagSuc > 0
        flagSuc = flagSuc -1;
        continue;
    end
    
    contendNs = find((regN(:,3)>0)&(regN(:,6)==0));
    
    if isempty(contendNs)
        continue;
    end
    
    % randomly pick
    sucN = contendNs(randi(length(contendNs)));
    
    % failed nodes register setup
    for i =1:length(contendNs)
        if contendNs(i) ~= sucN
            failN = contendNs(i);
            regN(failN, 6) = delayFail;
            regSuc(failN, regN(failN,2)+1) = regSuc(failN, regN(failN,2)+1)+1;
        end
    end
    
    % success node
    regN(sucN, 2) = regN(sucN, 2)+1;
    regN(sucN, 1) = (regN(sucN, 1) * (regN(sucN, 2)-1) + (t-regN(sucN, 4)))/regN(sucN, 2);
    regN(sucN, 3) = regN(sucN, 3) -1;
    regN(sucN, 4) = 0;
    
    regSuc(sucN, regN(sucN,2)) = regSuc(sucN, regN(sucN,2)) + 2;
    flagSuc = delaySuc;
end

%         numNResHard = numN - length(find(Nloc(:,3)>0));
%         numNResSoft = numN - length(find(Nloc(:,3)>0))/2;
%         for n = 1:numN
%             if(Nloc(n,3)>0)
%                 resHard(ii, nn) = resHard(ii, nn) + 1/length(find(Nloc(:,3)==Nloc(n,3)));
%                 resSoft(ii, nn) = resSoft(ii, nn) + 2/length(find(Nloc(:,3)==Nloc(n,3)));
%             else
%                 resHard(ii, nn) = resHard(ii, nn) + 1/numNResHard;
%                 resSoft(ii, nn) = resSoft(ii, nn) + 1/numNResSoft;
%
%             end
%             resNoHet(ii, nn) = resNoHet(ii, nn) + 1/numN;
%         end
%         resHard(ii,nn) = resHard(ii,nn) / numN;
%         resSoft(ii,nn) = resSoft(ii,nn) / numN;
%         resNoHet(ii,nn) = resNoHet(ii,nn) / numN;

toc

figure(1)
final_resHard = mean(resHard);
final_resSoft = mean(resSoft);
final_resNoHet = mean(resNoHet);
plot(final_resHard)
hold on;
plot(final_resSoft)
hold on;
plot(final_resNoHet)

%% results
%
% * ITEM1
% * ITEM2
%
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