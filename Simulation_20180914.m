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

tUnit = 1000;
numRes = maxT / tUnit;
Res = zeros(2,numRes);
trThres=0.1;

% for later use
nn = 1;

%% simulation start
numN = numNArray(nn);
numNtr = numN/2;
% Register Construction, 1:average Delay, 2: # delivered pkts, 3:#
% packets in a queue, 4:time of first packet, 5: next packet timer, 6: next trial timer
% Register Success: 1st #trials, 4:2nd...
regN = zeros(numN,6);
% 중간에 0으로 나누는 거 방지 위해 1부터 시작
regN(:,2) = 1;
GR = ceil(random('exp', ones(numN,1)*pktGenRate));
%regSuc = (-1)*ones(numN, maxT/min(GR));

t = 0;
flagSuc =0;
while (t<maxT)
    
    t= t+1;
    if flagSuc > 0
        regN(:,5) = regN(:,5)+ (regN(:,5)==0).* GR;
        regN(:,3) = regN(:,3) + (regN(:,5)<flagSuc);
        regN(:,4) = regN(:,4) + t*((regN(:,5)<flagSuc)&(regN(:,3)==1));
    
        regN(:,5) = (regN(:,5)-flagSuc);
        regN(:,5) = (regN(:,5)>0).* regN(:,5);
    
        regN(:,6) = (regN(:,6)-flagSuc);
        regN(:,6) = (regN(:,6)>0).* regN(:,6);
        t = t+flagSuc -1;
        flagSuc=0;
        continue;
    end
    
    
    
    % packet generation
    regN(:,5) = regN(:,5)+ (regN(:,5)==0).* GR;
    regN(:,3) = regN(:,3) + (regN(:,5)==1);
    regN(:,4) = regN(:,4) + t*((regN(:,5)==1)&(regN(:,3)==1));
    regN(:,5) = (regN(:,5)>0).* (regN(:,5)-1);
    regN(:,6) = (regN(:,6)>0).* (regN(:,6)-1);
    
    % WiFi Transition
    trMask = (regN(1:numNtr,3)>0)&(regN(1:numNtr,6)==0)&(rand(numNtr,1)<trThres);
    
    regN(1:numNtr, 2) = regN(1:numNtr, 2) + trMask;
    regN(1:numNtr, 1) = (regN(1:numNtr, 1) .* (regN(1:numNtr, 2)-1) + (trMask==0).*regN(1:numNtr, 1)+ (trMask==1).*(t+2000-regN(1:numNtr, 4)))./regN(1:numNtr, 2);
    regN(1:numNtr, 3) = regN(1:numNtr, 3) - trMask;
    regN(1:numNtr, 4) = (trMask==0).*regN(1:numNtr, 4);
    
    
    
    % success timer setup
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
 %           regSuc(failN, regN(failN,2)+1) = regSuc(failN, regN(failN,2)+1)+1;
        end
    end
    
    % success node
    regN(sucN, 2) = regN(sucN, 2)+1;
    regN(sucN, 1) = (regN(sucN, 1) * (regN(sucN, 2)-1) + (t-regN(sucN, 4)))/regN(sucN, 2);
    regN(sucN, 3) = regN(sucN, 3) -1;
    regN(sucN, 4) = 0;
    
  %  regSuc(sucN, regN(sucN,2)) = regSuc(sucN, regN(sucN,2)) + 2;
    flagSuc = delaySuc;
    
    Res(1,fix(t/tUnit)+1) = mean(regN(1:numN/2,1));
    Res(2,fix(t/tUnit)+1) = mean(regN(numN/2+1:numN,1));
end
toc

figure(1);
plot(Res(1,:),'-.b');
hold on;
plot(Res(2,:), '-.r');

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