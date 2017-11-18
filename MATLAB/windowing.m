function [BayesResult,Wforce] = windowing(data,FORCE,WindowSize_ms,fs,overlap, bins, alpha, beta, MVC)
% data [channels x samples]

% ////// TEST DATA ///// DELETE //////
% data = Test_data(:,1:22000);
% FORCE = Cal_force;
% fs = 1200;
% WindowSize_ms = 0.200;
% overlap = 70;
% /////////////////////////////////////

% -- Test mode -- %
buffering = 1; 

% -- Bayesian parameters -- %
% bins = 200;
% alpha = 10^-20;
% beta = 10^-500;
% MVC = 0.8;

% -- Windowing operation --%
WindowSize = floor((WindowSize_ms * fs));
newWindowSize = WindowSize;
factor = 1;
while (newWindowSize > 30)
    factor = factor +1;
    newWindowSize = floor(WindowSize / factor);
end

WindowIncrement = WindowSize - round(WindowSize*overlap/100);
LengthData = size(data,2);
LENGTH_FEAT = floor((LengthData - WindowSize )/WindowIncrement + 1);
CurrIx = 1;
BAY = [];
Wbayes = [];
Tbayes = zeros(size(data,1),factor);
temp_data = zeros(size(data,1),WindowSize);
indx = 1;
BayesianSignal = []; %zeros(size(data,1),floor(LENGTH_FEAT/factor));
[pri,sig,param] = loadbayes(fs,bins,alpha,beta,MVC); %reset pri

for p = 1:size(data,1)
    post(:,p) = pri;
end

while (CurrIx+WindowSize) <= LengthData
    
    CurrWindow = CurrIx:WindowSize+CurrIx-1;
    
    
    % --------------- Where magic happens ----------------------- %
    if (buffering == 1)
        for i = 1:size(data,1) % electrodes      
            temp_data(i,:) = data(i,CurrWindow)/param.sigmaMVC;
            % Second buffer to do not overload 'BayesFilter' requirements
            % (samples<40) and to reduce CPU time
            Wbayes = [];
            for k = 1:factor
                buffer = temp_data(i,indx:(floor(size(temp_data,2)/factor))+indx-1);
                [Bufbayes, post(:,i)] = BayesFilter(buffer, post(:,i), sig, param);
                Wbayes = [Wbayes Bufbayes]; % Wbayes [channels x factor]
                indx = indx+floor(size(temp_data,2)/factor);
            end
            indx = 1;
            Tbayes(i,:) = Wbayes;
        end
        BayesianSignal = [BayesianSignal Tbayes];
        
    elseif buffering == 0
         for i = 1:size(data,1) % electrodes 
            temp_data(i,:) = (data(i,CurrWindow)) / param.sigmaMVC;
            for j = 1:size(temp_data,2)
                [Wbayes, pri] = BayesFilter(temp_data(i,j), pri, sig, param);
            end
            Tbayes(i,:) = Wbayes;
        end
        BayesianSignal = [BayesianSignal Tbayes];
    
    elseif buffering == 2
         temp_data(1,:) = data(1,CurrWindow)/param.sigmaMVC;
            % Second buffer to do not overload 'BayesFilter' requirements
            % (samples<40) and to reduce CPU time
            Wbayes = [];
            for k = 1:factor
                buffer = temp_data(1,indx:(floor(size(temp_data,2)/factor))+indx-1);
                [Bufbayes, pri] = BayesFilter(buffer, pri, sig, param);
                Wbayes = [Wbayes Bufbayes]; % Wbayes [channels x factor]
                indx = indx+floor(size(temp_data,2)/factor);
            end
            indx = 1;
            Tbayes(1,:) = Wbayes;
            BayesianSignal = [BayesianSignal Tbayes];
    end    
    
    % ---------------------------------------------------------- %
    
    CurrIx = CurrIx+WindowIncrement;
end

if buffering == 2
BayesResult = BayesianSignal(1,:);
else 
BayesResult = BayesianSignal;
end

FOR = double(FORCE);
Wforce = ricampiona(FOR(1,:),size(BayesResult,2),'linear');
%plot(Wforce)
% hold on;plot(BayesResult(2,:),'.', 'linewidth',1)
% title('Bayesian signal');
disp(['windowing done']);
% hold on
% plot(Wforce)
%figure;plot(Cal_data')
%close figure;plot(bayesSTD')
end