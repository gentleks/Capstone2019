close all;
clear variables;

tic
%%%%%%%%%%%%%%%%
% Parameters   %
%%%%%%%%%%%%%%%%

iter = 1;

maxX = 100;
TransRan = 10;

numNStart = 10000;
numNStep = 100;
numNMax = 10000;
numNArray = numNStart:numNStep:numNMax;
NN = length(numNArray);
nn=1;

% delay parameter
delayFail = 10;
delaySuc = 200;
delayWiFi = 2000;

%simulation time
maxT = 1000*60*60*100;
pktGenRate = 1000*60*60;

tUnit = 1000;
numRes = maxT / tUnit;
Res = zeros(3,numRes);
trThres=0.1;


%% simulation start
numN = numNArray(nn);
numNtr = numN/2;
% Register Construction, 1:average Delay, 2: # delivered pkts, 3:#
% packets in a queue, 4:time of first packet, 5: next packet timer, 6: next trial timer
% Register Success: 1:decision(Stay/WiFi), 2:average, 3:current trials


for Type = 1:3
    
    regN = zeros(numN,5);
    
    GR = ceil(random('exp', ones(numN,1)*pktGenRate));
    regSuc = zeros(numN,3);
    
    t = 0;
    flagSuc =0;
    
    
    while (t<maxT)
        
        t= t+1;
        if flagSuc > 0
            %성공하면 패킷 generation 해주고, 나머지 패킷들의 timer setup
            regN(:,5) = regN(:,5)+ (regN(:,5)==0).* GR;
            regN(:,3) = regN(:,3) + (regN(:,5)<flagSuc);
            regN(:,4) = regN(:,4) + t*((regN(:,5)<flagSuc)&(regN(:,3)==1));
            
            regN(:,5) = max(0,(regN(:,5)-flagSuc));
            t = t+flagSuc -1;
            flagSuc=0;
            continue;
        end
        
        % packet generation
        regN(:,5) = regN(:,5)+ (regN(:,5)==0).* GR;
        regN(:,3) = regN(:,3) + (regN(:,5)==1);
        regN(:,4) = regN(:,4) + t*((regN(:,5)==1)&(regN(:,3)==1));
        
        % WiFi Transition....
        % trMask = (regN(1:numNtr,3)>0)&(rand(numNtr,1)<trThres);
        switch (Type)
            case 1
                for n = 1:numNtr
                    if regSuc(n,1) ==1 && regN(n,3)>0
                        regN(n, 2) = regN(n,2) +1;
                        regN(n, 1) = (regN(n,1) * (regN(n,2)-1) + delayWiFi)/ regN(n, 2);
                        regN(n, 3) = regN(n,3) -1;
                        regN(n, 4) = 0;
                    end
                end
            case 2
                trMask = regN(1:numNtr,3)>0;
                
                regN(1:numNtr, 2) = regN(1:numNtr, 2) + trMask;
                regN(1:numNtr, 1) = (regN(1:numNtr, 1) .* (regN(1:numNtr, 2)-1) + (trMask==0).*regN(1:numNtr, 1)+ (trMask==1).*delayWiFi)./ max(1,regN(1:numNtr, 2));
                regN(1:numNtr, 3) = regN(1:numNtr, 3) - trMask;
                regN(1:numNtr, 4) = (trMask==0).*regN(1:numNtr, 4);
            otherwise
        end
        
        regN(:,5) = max(0, regN(:,5)-1);
        
        % success timer setup
        contendNs = find((regSuc(:,1)==0) & (regN(:,3)>0));
        
        if isempty(contendNs)
            continue;
        end

        % randomly pick
        sucN = contendNs(randi(length(contendNs)));
        
        %fail register setup
        if Type ==1
            for i =1:length(contendNs)
                if contendNs(i) ~= sucN && regSuc(contendNs(i),1)==0
                    regSuc(contendNs(i), 3) = regSuc(contendNs(i), 3)+1;
                end
            end
        end
        
        % success node
        regN(sucN, 2) = regN(sucN, 2)+1;
        regN(sucN, 1) = (regN(sucN, 1) * (regN(sucN, 2)-1) + (t-regN(sucN, 4)))/regN(sucN, 2);
        regN(sucN, 3) = regN(sucN, 3) -1;
        regN(sucN, 4) = 0;
        
        if Type ==1
            if regSuc(sucN,1) == 1
                regSuc(sucN,3) = regSuc(sucN,2);
            end
            regSuc(sucN, 2) = (regSuc(sucN, 2) * (regN(sucN, 2)-1) + regSuc(sucN,3) )/regN(sucN, 2);
            regSuc(sucN, 3) = 0;
            if regSuc(sucN,2) * delaySuc > delayWiFi && sucN <=numNtr
                regSuc(sucN,1) = 1;
            end
        end
        flagSuc = delaySuc;
        
        Res(Type,fix(t/tUnit)+1) = sum(regN(1:numNtr,1).*regN(1:numNtr,2))/sum(regN(1:numNtr,2));
    end
end
toc

figure(1);
plot(Res(1,:),'-.r');
hold on;
plot(Res(2,:), '-.g');
hold on;
plot(Res(3,:), '-.b');

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