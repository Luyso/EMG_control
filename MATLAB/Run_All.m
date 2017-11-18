% 
% mean(Cortafuegos) % Are you sure you want to run all?
% %% Clean data
% addpath('Data', 'Michele');
% Subject = 6;
% [emg,force] = loadsubject(Subject);

% plot_array(activation')
% figure
% plot_array(force')

%%
% --- Bayesian filter --%

% Choose Subject
addpath('Data','Michele');
clear D1 D2 D3 D23 RMSE1 RMSE2 RMSE3 RMSE23 mape1 mape2 mape3 mape23
% Load data
StartSubjects = 1;
TotalSubjects = 1;
MAXDOF = 4;
%k = 6; % Number of synergies
MAXK = 2;
MINK = 2;
electrodes = 10;
samplingfreq = 2000;
% bins = 180;
% alpha = 10^-40;
% beta = 10^-150;
% MVC = 0.8;
disp(['Estimated Time = ',num2str(((TotalSubjects+1)*0.5)*(MAXK-MINK+1)), ' mins'])

% SubBayes = struct('K',[],'INFO',...
%     {{electrodes,alpha,beta,MAXDOF}});

GAIN = 5000; % Maximumm amplitud ~1.5. Cambiar la ganancia hace
% que la señal positiva y negativa se invierta para algún DOF

for Subject = StartSubjects:TotalSubjects
    [emg,force] = loadsubject(Subject);
    for k = MINK:MAXK
        disp(['Subject ',num2str(Subject), ' - k = ',num2str(k)]);
        Synergy = zeros(k*MAXDOF,electrodes);
        % DOF #1 (Index)    {0 - 50001}         3 repetitions
        % DOF #2 (Middle)   {100000 - 150000}   3 repetitions
        % DOF #3 (Ring)     {200000 - 250000}   3 repetitions
        % DOF #4 (little)   {295000 - 345000}   3 repetitions
        
        for DOF = 1:MAXDOF
            
            if DOF == 1
                start = 1;
                finish = 50001;
            elseif DOF == 2
                start = 100000;
                finish = 150000;
            elseif DOF == 3
                start = 200000;
                finish = 250000;
            elseif DOF == 4
                start = 295000;
                finish = 345000;
            end
            
            samples = finish - start;
            signal = zeros(samples,electrodes);
            
            % load new data
            for i = 1:electrodes
                for j = 1:samples
                    signal(j,i) = emg(j+start,i);
                end
            end
            signal = signal * GAIN;
            
            % -- Windowing parameters -- %
            samplingfreq = 1200 ;
            wTime        = 0.15 ; % Windows time in seconds
            over         = 50   ; % Windows overlap in percentage.
            
            bins = 100;
            alpha = 10^-30;
            beta = 10^-280;
            MVC = 0.8;
            % -- Perform bayessian filtering -- %
            clear bayesSTD
            disp(['Bayes filter of DOF #', num2str(DOF)]);
            [bayesSTD force_res(:,DOF)] = windowing(signal',force(start:finish,DOF)',wTime,samplingfreq,over,bins,alpha,beta,MVC);
        

            %
            %           figure
%            plot(bayesSTD') 
            %
            % ------------- NMF algorithm -------------
            
            disp(['Estimate Synergies of DOF #', num2str(DOF)]);
            % Start calibration
            Z = bayesSTD; % Each column of Z is a sample vector
            
            if max(max(isnan(Z))) == 1
                Z(isnan(Z)) = 0.1;
                NaN = 1
                disp(['There are some NaN in DOF #', num2str(DOF), '!!!'])
            end
            rank(Z);
            if k==size(Z,1)
                W = ones(k,k);
                H = Z;
                fprintf('#SYN==#CHAN\n\t W = ones(%d,%d);\n\t H = Z;\n',k,k) 
            else
                [W,H] = nnmf(Z,k,'replicates',10);
            end
            if DOF == 1
                ControlDOF1 = zeros(size(H));
                BayesDOF1 = bayesSTD;
                for j = 1:k
                    Synergy(j,:) = W(:,j);
                    ControlDOF1(j,:) = H(j,:);
                end
                
            elseif DOF == 2
                ControlDOF2 = zeros(size(H));
                BayesDOF2 = bayesSTD;
                for j = 1:k
                    Synergy(j+k,:) = W(:,j);
                    ControlDOF2(j,:) = H(j,:);
                end
                
            elseif DOF == 3
                ControlDOF3 = zeros(size(H));
                BayesDOF3 = bayesSTD;
                for j = 1:k
                    Synergy(j+(2*k),:) = W(:,j);
                    ControlDOF3(j,:) = H(j,:);
                end
                
            elseif DOF == 4
                ControlDOF4 = zeros(size(H));
                BayesDOF4 = bayesSTD;
                for j = 1:k
                    Synergy(j+(3*k),:) = W(:,j);
                    ControlDOF4(j,:) = H(j,:);
                end
                
            end
            disp('Done!')
            % ------------- NMF algorithm END -------------
            
            
        end
        SubBayes(Subject).K(k).SYN = Synergy;
        SubBayes(Subject).K(k).CAL = [ControlDOF1;ControlDOF2;ControlDOF3];
        if MAXDOF == 4
            SubBayes(Subject).K(k).CAL = [ControlDOF1;ControlDOF2;ControlDOF3;ControlDOF4];
        end
        
        disp('Calibration complete.')
        
        % ------------- Regressor -------------
        maxRang = [-1 15];
        disp('Calculate regressor parameters...')
        % DOF #1
        yd = force_res(:,1);
        x_tr = ControlDOF1(1,:)';
        for i = 2:k
            x_tr = [x_tr ControlDOF1(i,:)'];
        end
        X = [zeros(size(x_tr,1),1) x_tr];
        beta1_1 = regress(yd,X);
        Csum = beta1_1(1);
        for j = 1:k
            Ctot(j,:) = beta1_1(j+1) * ControlDOF1(j,:);
            Csum = Csum + Ctot(j,:);
        end
        figure;subplot(1,MAXDOF,1);hold on;plot(yd);plot(Csum); title('DOF #1')
        axis([1 size(bayesSTD,2) maxRang])
        %     x_tr = [ControlDOF1(1,:)' ControlDOF1(2,:)'];
        %     X = [ones(size(x_tr,1),1) x_tr];
        %     beta1_1 = regress(yd,X);
        % DOF #2
        yd = force_res(:,2);
        x_tr = ControlDOF2(1,:)';
        for i = 2:k
            x_tr = [x_tr ControlDOF2(i,:)'];
        end
        X = [zeros(size(x_tr,1),1) x_tr];
        beta2_1 = regress(yd,X);
        Csum = beta2_1(1);
        for j = 1:k
            Ctot(j,:) = beta2_1(j+1) * ControlDOF2(j,:);
            Csum = Csum + Ctot(j,:);
        end
            subplot(1,MAXDOF,2);hold on;plot(yd);plot(Csum);title('DOF #2')
        axis([1 size(bayesSTD,2) maxRang])
        %DOF #3
        yd = force_res(:,3);
        x_tr = ControlDOF3(1,:)';
        for i = 2:k
            x_tr = [x_tr ControlDOF3(i,:)'];
        end
        X = [zeros(size(x_tr,1),1) x_tr];
        beta3_1 = regress(yd,X);
        
        Csum = beta3_1(1);
        for j = 1:k
            Ctot(j,:) = beta3_1(j+1) * ControlDOF3(j,:);
            Csum = Csum + Ctot(j,:);
        end
            subplot(1,MAXDOF,3);hold on;plot(yd);plot(Csum);title('DOF #3')
        axis([1 size(bayesSTD,2) maxRang])
        
        % DOF #4
        if MAXDOF == 4
            yd = force_res(:,4);
            x_tr = ControlDOF4(1,:)';
            for i = 2:k
                x_tr = [x_tr ControlDOF4(i,:)'];
            end
            X = [zeros(size(x_tr,1),1) x_tr];
            beta4_1 = regress(yd,X);
        
        Csum = beta4_1(1);
        for j = 1:k
            Ctot(j,:) = beta4_1(j+1) * ControlDOF4(j,:);
            Csum = Csum + Ctot(j,:);
        end
            subplot(1,MAXDOF,4);hold on;plot(yd);plot(Csum);title('DOF #4')
        axis([1 size(bayesSTD,2) maxRang])
        end
        disp('Done!')
        %close all
        % ------------- Regressor END ------------- %
        
        
        % -------------------- Calibration END ----------------------------- %
        % --------------------- OFFLINE TEST ------------------------------- %
        disp(['Offline Test start for Subject ', num2str(Subject)])
        
        % TDOF #1 (Index) {50000 - 100000} 3 reps
        % TDOF #2 (Middle) {165000 - 210000} 3 reps
        % TDOF #3 (Ring) {246000 - 296000} 3 reps
        % TDOF #4 (Middle+Ring) {735000 - 780000} 3 reps
        
        TESTDOF = 6;
        for TDOF = 1:TESTDOF
            
            if TDOF == 1
                Tstart = 50000;
                Tfinish = 100000;
            elseif TDOF == 2
                Tstart = 150000;
                Tfinish = 200000;
            elseif TDOF == 3
                Tstart = 250000;
                Tfinish = 300000;
            elseif TDOF == 4 %DOF 4
                Tstart = 345000;
                Tfinish = 395000;
            elseif TDOF == 5 % DOF 1 and 4
                Tstart = 600000;
                Tfinish = 650000;
            elseif TDOF == 6 % DOF 2 and 3
                Tstart = 730000;
                Tfinish = 780000;
            end
            
            Tsamples = Tfinish - Tstart;
            Tsignal = zeros(Tsamples,electrodes);
            
            % load new data
            for i = 1:electrodes
                for j = 1:Tsamples
                    Tsignal(j,i) = emg(j+Tstart,i);
                end
            end
            Tsignal = Tsignal * GAIN; % Maximumm amplitud ~1.5
            %
            %figure
            %plot (Tsignal)
            % -- Windowing parameters -- %
            % -- Perform bayessian filtering -- %
            clear TbayesSTD
            disp(['Bayes filter of TDOF #', num2str(TDOF)]);
            [TbayesSTD Tforce_res(:,DOF)] = windowing(Tsignal',force(Tstart:Tfinish,DOF)',wTime,samplingfreq,over,bins,alpha,beta,MVC);
            
            disp('Done!');
            %close all
            
            %%%%%%%%%%%%%%%%%%%%%%% END Bayessian Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% NMF algorithm %%
            
            disp('Estimating Control signals...');
            TZ = TbayesSTD; % Each column of Z is a sample vector
            Control = pinv(Synergy)' * TZ; % Using Moore-Penros alg.
%             Control =  Synergy * pinv(TZ)';            
%             Control = pinv(Control)';
            
            disp('Done!');
            
            % Reconstruct the signal
            for i = 1:k*MAXDOF
                C(i,:) = Control(i,:);
            end
            ca1 = C(1,:)*beta1_1(2);
            ca2 = C(k+1,:)*beta2_1(2);
            ca3 = C(1+2*k,:)*beta3_1(2);
            for i = 2:k
                ca1 = ca1+C(i,:)*beta1_1(i+1);
            end
            for i = 2:k
                ca2 = ca2+C(i+k,:)*beta2_1(i+1);
            end
            for i = 2:k
                ca3 = ca3+C(i+k*2,:)*beta3_1(i+1);
            end
            Sum1 = beta1_1(1) + ca1;
            Sum2 = beta2_1(1) + ca2;
            Sum3 = beta3_1(1) + ca3;
            Sum = [Sum1;Sum2;Sum3];
            % -- Resample forces -- %
            FOR = double(force);
            clear Tforce_res
            for i=1:MAXDOF
                Tforce_res(:,i) = ricampiona(FOR(Tstart:Tfinish,i),size(Sum,2),'linear');
            end
            
            For = [Tforce_res(:,1);Tforce_res(:,2);Tforce_res(:,3)];
            
            if MAXDOF == 4
                ca4 = C(1+3*k,:)*beta4_1(2);
                for i = 2:k
                    ca4 = ca4+C(i+3*k,:)*beta4_1(i+1);
                end
                Sum4 = beta4_1(1) + ca4;
                Sum = [Sum1;Sum2;Sum3;Sum4];
                For = [Tforce_res(:,1);Tforce_res(:,2);...
                    Tforce_res(:,3);Tforce_res(:,4)];
            end
            SubBayes(Subject).K(k).TDOF(TDOF).RawCONT = Control;
            SubBayes(Subject).K(k).TDOF(TDOF).CONT = Sum;
            SubBayes(Subject).K(k).TDOF(TDOF).BAY = TbayesSTD;
            MVC_Force = max(max(For));
            Tol = 0.05; % Increasing the tolerance may decrease the euclidian distance.
            Err = Tol*MVC_Force; %
            disp(['Calculate errors for TDOF #', num2str(TDOF)]);
            for i= 1:MAXDOF
                var1 = Sum(i,:);
                var2 = Tforce_res(:,i)';
                %                 var1(var1<0) = Err;
                %                 var2(var2<0) = Err;
                
                if TDOF == 1
                    RMSE1(Subject,i) = sqrt(mean((var1-var2).^2)); % Root mean squared error
                    mape1(Subject,i) = mean((var2+var1)/var2); % mean absolute percentage error
                    D1(Subject,i) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2)); % euclidian distance
                    RS1(Subject,i) = dAvella(var1,var2);
                    SubBayes(Subject).K(k).ERR.D(TDOF,i) = D1(Subject,i);
                    SubBayes(Subject).K(k).ERR.RMSE(TDOF,i) = RMSE1(Subject,i);
                    SubBayes(Subject).K(k).ERR.RS(TDOF,i) = RS1(Subject,i);
                elseif TDOF == 2
                    RMSE2(Subject,i) = sqrt(mean((var1-var2).^2)); % Root mean squared error
                    mape2(Subject,i) = mean((var2+var1)/var2); % mean absolute percentage error
                    D2(Subject,i) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2)); % euclidian distance
                    RS2(Subject,i) = dAvella(var1,var2);
                    SubBayes(Subject).K(k).ERR.D(TDOF,i) = D2(Subject,i);
                    SubBayes(Subject).K(k).ERR.RMSE(TDOF,i) = RMSE2(Subject,i);
                    SubBayes(Subject).K(k).ERR.RS(TDOF,i) = RS2(Subject,i);
                elseif TDOF == 3
                    RMSE3(Subject,i) = sqrt(mean((var1-var2).^2)); % Root mean squared error
                    mape3(Subject,i) = mean((var2+var1)/var2); % mean absolute percentage error
                    D3(Subject,i) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2)); % euclidian distance
                    RS3(Subject,i) = dAvella(var1,var2);
                    SubBayes(Subject).K(k).ERR.D(TDOF,i) = D3(Subject,i);
                    SubBayes(Subject).K(k).ERR.RMSE(TDOF,i) = RMSE3(Subject,i);
                    SubBayes(Subject).K(k).ERR.RS(TDOF,i) = RS3(Subject,i);
                elseif TDOF == 4
                    RMSE4(Subject,i) = sqrt(mean((var1-var2).^2)); % Root mean squared error
                    mape4(Subject,i) = mean((var2+var1)/var2); % mean absolute percentage error
                    D4(Subject,i) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2)); % euclidian distance
                    RS4(Subject,i) = dAvella(var1,var2);
                    SubBayes(Subject).K(k).ERR.D(TDOF,i) = D4(Subject,i);
                    SubBayes(Subject).K(k).ERR.RMSE(TDOF,i) = RMSE4(Subject,i);
                    SubBayes(Subject).K(k).ERR.RS(TDOF,i) = RS4(Subject,i);
                elseif TDOF == 6
                    RMSE23(Subject,i) = sqrt(mean((var1-var2).^2)); % Root mean squared error
                    mape23(Subject,i) = mean((var2+var1)/var2); % mean absolute percentage error
                    D23(Subject,i) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2)); % euclidian distance
                    RS23(Subject,i) = dAvella(var1,var2);
                    SubBayes(Subject).K(k).ERR.D(TDOF,i) = D23(Subject,i);
                    SubBayes(Subject).K(k).ERR.RMSE(TDOF,i) = RMSE23(Subject,i);
                    SubBayes(Subject).K(k).ERR.RS(TDOF,i) = RS23(Subject,i);
                elseif TDOF == 5
                    RMSE14(Subject,i) = sqrt(mean((var1-var2).^2)); % Root mean squared error
                    mape14(Subject,i) = mean((var2+var1)/var2); % mean absolute percentage error
                    D14(Subject,i) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2)); % euclidian distance
                    RS14(Subject,i) = dAvella(var1,var2);
                    SubBayes(Subject).K(k).ERR.D(TDOF,i) = D14(Subject,i);
                    SubBayes(Subject).K(k).ERR.RMSE(TDOF,i) = RMSE14(Subject,i);
                    SubBayes(Subject).K(k).ERR.RS(TDOF,i) = RS14(Subject,i);
                end
                
            end
        end
    end
    disp(['Estimated Time of acomplishment = ',num2str(((TotalSubjects - Subject)*0.5)*(MAXK-MINK+1)), ' mins'])
    
end
disp('Ouuuu yeah!')


%% Plot Results of specific Subject

%% Choose Subject and Test DOF
% ----------------
Subject     =    1
TDOF        =    1
k           =    2
Classifier  =    0
% -----------------
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% /////////// Classifier test ////////////////////
SubBayes = Subjects40_4DOFS;
Bastards = [];
TotalError = zeros(MAXK,MAXDOF);
% for k = 1:6
%     disp([' Test of k = ',num2str(k)])
% for TDOF = 1:MAXDOF
    Match = 0;
    fail = 0;
    for Subject = 1:20
        
[emg,force] = loadsubject(Subject);
    if Subject == 8
        Subject = 9;
    end
if TDOF == 1
    Tstart = 50000;
    Tfinish = 100000;
elseif TDOF == 2
    Tstart = 150000;
    Tfinish = 200000;
elseif TDOF == 3
    Tstart = 250000;
    Tfinish = 300000;
elseif TDOF == 4
    Tstart = 345000;
    Tfinish = 395000;
elseif TDOF == 5 % DOF 1 and 4
    Tstart = 600000;
    Tfinish = 650000;
elseif TDOF == 6 % DOF 2 and 3
    Tstart = 730000;
    Tfinish = 780000;
end


% Reconstruct signal
Sum1 = SubBayes(Subject).K(k).TDOF(TDOF).CONT(1,:);
Sum2 = SubBayes(Subject).K(k).TDOF(TDOF).CONT(2,:);
Sum3 = SubBayes(Subject).K(k).TDOF(TDOF).CONT(3,:);
if MAXDOF == 4
    Sum4 = SubBayes(Subject).K(k).TDOF(TDOF).CONT(4,:);
end
FOR = double(force);
clear Tforce_res
for i=1:MAXDOF
    Tforce_res(:,i) = ricampiona(FOR(Tstart:Tfinish,i),size(Sum1,2),'linear');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////// Classifier ///////////////////
Syn = SubBayes(Subject).K(k).SYN;
Bayess = SubBayes(Subject).K(k).TDOF(TDOF).BAY;
RawCont = SubBayes(Subject).K(k).TDOF(TDOF).RawCONT;
error = [];
Active = 0;
if Classifier == 1
    data = 1; % data = 1 for ninapro
    for s = 1:MAXDOF
        [Terror,Fing] = classifier(Bayess,RawCont,Syn,s,0,k,data);
% // Treshold based classification // %
%          if Terror < 15
%             disp(['DOF ',num2str(s),': Active'])
%             Active = [Active s];
%          end
         error = [error Terror];
    end
    error
    [~,Active] = min(error);
    
end
    if Active == TDOF
        Match = Match+1;
    else
       %disp(['Fail!! of Subject ',num2str(Subject),' for DOF: ',num2str(TDOF),', K = ', num2str(k)])
       fail = [fail Subject];
    end

    end
    ClasErr(TDOF) = Match;
    
    if TDOF == 1
        Bastards1 = fail;
    elseif TDOF == 2
        Bastards2 = fail;
    elseif TDOF == 3
        Bastards3 = fail;
    end
% end
TotalError(k,:) = ClasErr';
% end
disp('DONE')

%% -- Plot Classifier Results -- %%
figure
bar(TotalError./40)
title('Classification Accuracy for single DOF moves')
xlabel('Number of synergies (k)')
ylabel('% Accuracy ')
legend('DOF 1','DOF 2','DOF 3','DOF 4')
axis ([0.5,6.5,0,1])
% ///////////////////////////////////////////////
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ///////////////////////////////////////////////

% -- Plot results -- %%

% ----------------
Subject     =    14
TDOF        =    5
k           =    2
% -----------------

figure('units', 'normalized', 'outerposition', [0 0.05 1 0.5])

for k = 1:6
[emg,force] = loadsubject(Subject);
    if Subject == 8
        Subject =9;
    end
if TDOF == 1
    Tstart = 50000;
    Tfinish = 100000;
elseif TDOF == 2
    Tstart = 150000;
    Tfinish = 200000;
elseif TDOF == 3
    Tstart = 250000;
    Tfinish = 300000;
elseif TDOF == 4
    Tstart = 345000;
    Tfinish = 395000;
elseif TDOF == 5 % DOF 1 and 4
    Tstart = 600000;
    Tfinish = 650000;
elseif TDOF == 6 % DOF 2 and 3
    Tstart = 730000;
    Tfinish = 780000;
end

% Reconstruct signal
Sum1 = SubBayes(Subject).K(k).TDOF(TDOF).CONT(1,:);
Sum2 = SubBayes(Subject).K(k).TDOF(TDOF).CONT(2,:);
Sum3 = SubBayes(Subject).K(k).TDOF(TDOF).CONT(3,:);
Sum = [Sum1;Sum2;Sum3];
if MAXDOF == 4
    Sum4 = SubBayes(Subject).K(k).TDOF(TDOF).CONT(4,:);
    Sum = [Sum1;Sum2;Sum3;Sum4];
end
FOR = double(force);
clear Tforce_res
for i=1:MAXDOF
    Tforce_res(:,i) = ricampiona(FOR(Tstart:Tfinish,i),size(Sum,2),'linear');
end
clear S
for i = 1:k*MAXDOF
    S(i,:) = SubBayes(Subject).K(k).SYN(i,:);
end
Range = [-5 max(max(max(Tforce_res),max(max(Sum))))+1];
lim = [0,11,1,20];
samples = size(Sum,2);

%DOF #1
subplot(MAXDOF,6,k);hold on;plot(Tforce_res(:,1));plot(Sum1, 'LineWidth',1.5);
axis([1 samples Range])
title([ 'K = ',num2str(k)])
if k == 1
ylabel('Index')
end
subplot(MAXDOF,6,k+6);hold on;plot(Tforce_res(:,2));plot(Sum2, 'LineWidth',1.5)
axis([1 samples Range])
if k+6 == 7
ylabel('Middle')
end
subplot(MAXDOF,6,k+(6*2));hold on;plot(Tforce_res(:,3));plot(Sum3, 'LineWidth',1.5)
axis([1 samples Range])

if k+(6*2) == 13
ylabel('Ring')
end
if MAXDOF == 4
subplot(MAXDOF,6,k+(6*3));hold on;plot(Tforce_res(:,4));plot(Sum4, 'LineWidth',1.5)
axis([1 samples Range])
if k+(6*3) == 19
ylabel('Little')
end
end
end
%% Plot errors
Subject = 37
close all
MINK = 1;
MAXK = 6;
clear TD1 TD2 TD3 TD4 TD5 TD6
for i = MINK:MAXK
    TD1(i,:) = SubBayes(Subject).K(i).ERR.RMSE(1,:);
    TD2(i,:) = SubBayes(Subject).K(i).ERR.RMSE(2,:);
    TD3(i,:) = SubBayes(Subject).K(i).ERR.RMSE(3,:);
    TD4(i,:) = SubBayes(Subject).K(i).ERR.RMSE(6,:);
    TD = max([TD1 TD2 TD3 TD4]);
    if MAXDOF == 4
    TD4(i,:) = SubBayes(Subject).K(i).ERR.RMSE(4,:);
    TD5(i,:) = SubBayes(Subject).K(i).ERR.RMSE(5,:);
    TD6(i,:) = SubBayes(Subject).K(i).ERR.RMSE(6,:);
    TD = max([TD1 TD2 TD3 TD4 TD5 TD6]);
    end
end

subplot(2,3,1)
bar(TD1./size(SubBayes(Subject).K(1).CAL,2))
title('Error for DOF #1 active')
xlabel('k')
subplot(2,3,2)
bar(TD2./size(SubBayes(Subject).K(1).CAL,2))
title('Error for DOF #2 active')
subplot(2,3,3)
bar(TD3./size(SubBayes(Subject).K(1).CAL,2))
title('Error for DOF #3 active')
subplot(2,3,4)
bar(TD4./size(SubBayes(Subject).K(1).CAL,2))
title('Error for DOFs #4 active')
if MAXDOF == 4
subplot(2,3,5)
bar(TD5./size(SubBayes(Subject).K(1).CAL,2))
title('Error for DOFs #1 and #4 active')
subplot(2,3,6)
bar(TD6./size(SubBayes(Subject).K(1).CAL,2))
title('Error for DOFs #2 and #3 active')
end
%bar(SubBayes(Subject).SYN(1,:));hold all;bar(SubBayes(2).SYN(1,:));
%% -- Plot error of active DOFs (1 to 4) -- %%
clear E1 E2 E3 E4 E5 E6 
for m = 1:TotalSubjects
    for n = 1:MAXDOF
    E1(m,n) = SubBayes(m).K(1).ERR.RS(n,n);
    end
end

for m = 1:TotalSubjects
    for n = 1:MAXDOF
    E2(m,n) = SubBayes(m).K(2).ERR.RS(n,n);
    end
end
for m = 1:TotalSubjects
    for n = 1:MAXDOF
    E3(m,n) = SubBayes(m).K(3).ERR.RS(n,n);
    end
end
for m = 1:TotalSubjects
    for n = 1:MAXDOF
    E4(m,n) = SubBayes(m).K(4).ERR.RS(n,n);
    end
end
for m = 1:TotalSubjects
    for n = 1:MAXDOF
    E5(m,n) = SubBayes(m).K(5).ERR.RS(n,n);
    end
end
for m = 1:TotalSubjects
    for n = 1:MAXDOF
    E6(m,n) = SubBayes(m).K(6).ERR.RS(n,n);
    end
end
% for i = 1:MAXDOF
% subplot(4,MAXDOF,i)
% bar(E1(:,i))
% title(['Error of DOF #', num2str(i),' when DOF #', num2str(i),' is active'])
% ylabel('k = 1')
% axis([1,TotalSubjects+1,0,10])
% subplot(4,MAXDOF,i+MAXDOF)
% bar(E2(:,i))
% ylabel('k = 2')
% axis([1,TotalSubjects+1,0,10])
% subplot(4,MAXDOF,i+MAXDOF*2)
% bar(E4(:,i))
% ylabel('k = 4')
% axis([1,TotalSubjects+1,0,10])
% subplot(4,MAXDOF,i+MAXDOF*3)
% bar(E6(:,i))
% ylabel('k = 6')
% axis([1,TotalSubjects+1,0,10])
% end

figure
E1(isnan(E1)) = 3;
E2(isnan(E2)) = 3;
E3(isnan(E3)) = 3;
E4(isnan(E4)) = 3;
E5(isnan(E5)) = 3;
E6(isnan(E6)) = 3;
means = [mean(E1);mean(E2);mean(E3);mean(E4);mean(E5);mean(E6)];
med = [median(E1);median(E2);median(E3);median(E4);median(E5);median(E6)];

bar(means)
title(' ERROR: Control estimated signal Vs force measured of only active DOFs.')
xlabel('Number of synergies (k)')
ylabel('mean error')
legend('DOF 1','DOF 2','DOF 3','DOF 4')

