clear force_res Ctot

% ////////////////////////////////////////////////////////////////////////////
% Run two different process: Calibration and Offline Test. 
% The total number of DOF's is set to 4 (index, middle, ring and little finger respectively)
% How to run:
% 1 - Choose the number of subjects (up to 40)
% 2 - Choose the number of synergies
% 3 - Choose the reescaling approach (reescaling.m)
% 4 - Run
% After the offline test is completed the program will stop. 
%
% The error results are save in the structure:
% NinaResults(reescaling).Subject(s).K(k).ERR.D(TASK,DOF) 
% NinaResults(reescaling).Subject(s).K(k).ERR.RMSE(TASK,DOF) 
% NinaResults(reescaling).Subject(s).K(k).ERR.RS(TASK,DOF) 
%
% DATASET:
% The data includes six tasks (flexions), six repetitions each: 
% 1- Index
% 2- Middle
% 3- Ring
% 4- Little
% 5- Index + Little
% 6- Ring + Little
% 
% The data is split into 3 repetitions for calibration and 3 repetitions for testing.
% ////////////////////////////////////////////////////////////////////////////



addpath('Data')
% Set NinaPro params
electrodes = 10;
samplingfreq = 2000;
option = 'Proportional';
gain = 1/5000;
plotear = 1;

if strcmp(option,'Proportional')
    scaleOpt = 1;
elseif strcmp(option,'Normalize')
    scaleOpt = 2;
elseif strcmp(option,'Rescale')
    scaleOpt = 3;
elseif strcmp(option,'MAX')
    scaleOpt = 4;
end


% Total number of DOF's
MAXDOF = 4;

% Chose number of subjects
StartSubject = 1;
FinalSubject = 1;

% Set number of synergies
MINK = 2;
MAXK = 2;

% Non-negative matrix factorization mode
% 'nmf' -> Regular NNMF
% 'snmf' -> Sparse NMF
% Algorithm = 'snmf';
% Sparsity = 0.3;


if (FinalSubject - StartSubject) > 1 || (MAXK - MINK) > 1
    plotear = 0;
end

disp(['Estimated Time = ',num2str(((FinalSubject+1)*0.5)*(MAXK-MINK+1)), ' mins'])

% ------ Calibration phase START ---------- %

for Subject = StartSubject:FinalSubject
    
    [emg,force] = loadsubject(Subject);
    
    for k = MINK:MAXK
        disp(['Subject ',num2str(Subject), ' - k = ',num2str(k)]);
        
        Synergy = zeros(k*MAXDOF,electrodes);
        

        
        for DOF = 1:MAXDOF
            
            [start, finish] = loadindexNINA(DOF,1);
            samples = finish - start;
            
            signal = emg(start:finish,1:electrodes);
            
            % Reescale emg signal before apply bayesian filtering

            [Cal_data,Cal_force] = reescaling(signal',force',DOF,[],0,option,gain);
            
            % -- Windowing parameters -- %
            samplingfreq = 1200 ;
            wTime        = 0.15 ; % Windows time in seconds
            over         = 50   ; % Windows overlap in percentage.
            
            % -- Bayes parameters -- %
            bins = 100;
            alpha = 10^-7;
            beta = 10^-250;
            MVC = 0.8;
            
            % -- Perform bayessian filtering -- %
            disp(['Bayes filter of DOF #', num2str(DOF)]);
            [bayesSTD force_res(:,DOF)] = windowing(Cal_data,force(start:finish,DOF)',wTime,samplingfreq,over,bins,alpha,beta,MVC);
            
            %% ------------- NMF algorithm START ------------- %%
            
            disp(['Estimate Synergies of DOF #', num2str(DOF)]);
            % Start calibration
            Z = bayesSTD; % Each column of Z is a sample vector
            
            % Check for NaN values
            if max(max(isnan(Z))) == 1
                Z(isnan(Z)) = 0.1;
                NaN = 1
                disp(['There are some NaN in DOF #', num2str(DOF), 'of subject  ', num2str(Subject),'!!!'])
            end
            
            opt = statset('MaxIter',10); %'Display','final');
            [W0,H0] = nnmf(Z,k,'replicates',5,'options',opt,'algorithm','mult');
            opt = statset('Maxiter',100);
            [W,H] = nnmf(Z,k,'w0',W0,'h0',H0,'options',opt,'algorithm','als');
            
            %                 if Algorithm == 'snmf'
            %                     options.beta = Sparsity;
            %                     options.eta = max(max(Z))^2;
            %                     [A,Y,numIter,tElapsed,finalResidual] = sparsenmfnnls(Z,k,options);
            %                     W = A;
            %                     H = Y;
            %                 end
            
            
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
        NinaResult(scaleOpt).Sub(Subject).K(k).SYN = Synergy;
        NinaResult(scaleOpt).Sub(Subject).K(k).CAL = [ControlDOF1;ControlDOF2;ControlDOF3];
        if MAXDOF == 4
            NinaResult(scaleOpt).Sub(Subject).K(k).CAL = [ControlDOF1;ControlDOF2;ControlDOF3;ControlDOF4];
        end
        
        disp('Calibration complete.')
        
        % ------------- Regression START ------------- %
        
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
        if plotear == 1
        figure;subplot(1,MAXDOF,1);hold on;plot(yd);plot(Csum); title('DOF #1')
        end
        %axis([1 size(bayesSTD,2) maxRang])
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
        if plotear == 1
        subplot(1,MAXDOF,2);hold on;plot(yd);plot(Csum);title('DOF #2')
        end
        %axis([1 size(bayesSTD,2) maxRang])
        
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
        if plotear == 1
        subplot(1,MAXDOF,3);hold on;plot(yd);plot(Csum);title('DOF #3')
        %axis([1 size(bayesSTD,2) maxRang])
        end
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
            if plotear == 1
            subplot(1,MAXDOF,4);hold on;plot(yd);plot(Csum);title('DOF #4')
            %axis([1 size(bayesSTD,2) maxRang])
            end
        end
        disp('Done!')
        %close all
        
        % ------------- Regressor END ------------- %
        
        
        % -------------------- Calibration phase END ----------------------------- %
        
        
        %% --------------------- OFFLINE TEST ------------------------------- %%
        disp(['Offline Test start for Subject ', num2str(Subject)])
        

        % TASK = 1;
        
        for TASK = 1:6
            
            [Tstart, Tfinish] = loadindexNINA(TASK,1);
            Tsamples = Tfinish - Tstart;
            Tsignal = emg(Tstart:Tfinish,1:electrodes);           
            [Test_data,Test_force] = reescaling(Tsignal',force',DOF,[],0,option,gain);
            %figure
            %plot (Test_data)
            
            % -- Perform Bayesian filtering -- %
            
            disp(['Bayes filter of TASK #', num2str(TASK)]);
            [TbayesSTD, ~] = windowing(Test_data,force(Tstart:Tfinish,DOF)',wTime,samplingfreq,over,bins,alpha,beta,MVC);
            disp('Done!');
            %close all
            
            %%%%%%%%%%%%%%%%%%%%%%% END Bayessian Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% NMF algorithm %%
            
            disp('Estimating Control signals...');
            TZ = TbayesSTD; % Each column of TZ is a sample vector
            Control = pinv(Synergy)' * TZ; % Using Moore-Penros alg.
            disp('Done!');
            
            % Reconstruct signal
            C = Control;
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
            
            NinaResult(scaleOpt).Sub(Subject).K(k).TASK(TASK).RawCONT = Control;
            NinaResult(scaleOpt).Sub(Subject).K(k).TASK(TASK).CONT = Sum;
            NinaResult(scaleOpt).Sub(Subject).K(k).TASK(TASK).BAY = TbayesSTD;
            
            % MVC_Force = max(max(For));
            % Tol = 0.05; % Increasing the tolerance may decrease the euclidian distance
            % Err = Tol*MVC_Force;
            
            disp(['Calculate errors for TASK #', num2str(TASK)]);
            for j = 1:4
                var1 = Sum(j,:);
                var2 = Tforce_res(:,j)';
                % var1(var1<0) = Err;
                % var2(var2<0) = Err;
                energy = (sum((var1).^2));
                NinaResult(scaleOpt).Sub(Subject).K(k).ERR.D(TASK,j) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2));
                NinaResult(scaleOpt).Sub(Subject).K(k).ERR.RMSE(TASK,j) = sqrt(sum((var1-var2).^2)/energy);
                NinaResult(scaleOpt).Sub(Subject).K(k).ERR.RS(TASK,j) = Rsquare(var1, var2);
            end
            disp('Done!');
        end
        disp(['OFFLINE test for subject ', num2str(Subject), ' complete!']);
    end
    disp(['ETA = ',num2str(((FinalSubject - Subject)*0.5)*(MAXK-MINK+1)), ' mins'])
end
disp('Ouuuu yeah!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ///////////////////////////////////////////////

%% Plot Results for a specific Subject and TASK

% ----------------
Subject     =    37
TASK        =    6
% k           =    2
% -----------------

figure('units', 'normalized', 'outerposition', [0 0.05 1 0.5])

for k = MINK:MAXK
    [emg,force] = loadsubject(Subject);
    if Subject == 8
        Subject =9;
    end
    [Tstart,Tfinish] = loadindexNINA(TASK,1);
    
    % Reconstruct signal
    Sum1 = NinaResult(scaleOpt).Sub(Subject).K(k).TASK(TASK).CONT(1,:);
    Sum2 = NinaResult(scaleOpt).Sub(Subject).K(k).TASK(TASK).CONT(2,:);
    Sum3 = NinaResult(scaleOpt).Sub(Subject).K(k).TASK(TASK).CONT(3,:);
    Sum = [Sum1;Sum2;Sum3];
    if MAXDOF == 4
        Sum4 = NinaResult(scaleOpt).Sub(Subject).K(k).TASK(TASK).CONT(4,:);
        Sum = [Sum1;Sum2;Sum3;Sum4];
    end
    FOR = double(force);
    clear Tforce_res
    for i=1:MAXDOF
        Tforce_res(:,i) = ricampiona(FOR(Tstart:Tfinish,i),size(Sum,2),'linear');
    end
    clear S
    for i = 1:k*MAXDOF
        S(i,:) = NinaResult(scaleOpt).Sub(Subject).K(k).SYN(i,:);
    end
    Range = [-5 max(max(max(Tforce_res),max(max(Sum))))+1];
    Range_Mid = [-5 max(max(max(Tforce_res),25))];
    lim = [0,11,1,20];
    samples = size(Sum,2);
    
    %DOF #1
    subplot(MAXDOF,MAXK,k);hold on;plot(Tforce_res(:,1), 'LineWidth',2.5);plot(Sum1, 'LineWidth',2.5);
    axis([1 samples Range])
    title([ 'K = ',num2str(k)])
    if k == 1
        ylabel('Index')
    end
    subplot(MAXDOF,MAXK,k+MAXK);hold on;plot(Tforce_res(:,2),'LineWidth',2.5);plot(Sum2, 'LineWidth',2.5)
    axis([1 samples Range_Mid])
    if k+6 == 7
        ylabel('Middle')
    end
    subplot(MAXDOF,MAXK,k+(MAXK*2));hold on;plot(Tforce_res(:,3),'LineWidth',2.5);plot(Sum3, 'LineWidth',2.5)
    axis([1 samples Range])
    
    if k+(6*2) == 13
        ylabel('Ring')
    end
    if MAXDOF == 4
        subplot(MAXDOF,MAXK,k+(MAXK*3));hold on;plot(Tforce_res(:,4),'LineWidth',2.5);plot(Sum4, 'LineWidth',2.5)
        axis([1 samples Range])
        if k+(6*3) == 19
            ylabel('Little')
        end
    end
end


