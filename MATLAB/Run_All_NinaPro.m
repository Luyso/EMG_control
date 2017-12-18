clear C D1 D2 D3 D23 RMSE1 RMSE2 RMSE3 RMSE23 mape1 mape2 mape3 mape23 force_res Ctot Tforce_res

% Set NinaPro params
electrodes = 10;
samplingfreq = 2000;
GAIN = 5000; % Set

% Total number of DOF's
MAXDOF = 4;

% Chose number of subjects
StartSubject = 37;
FinalSubject = 37;

% Non-negative matrix factorization mode
% 'nmf' -> Regular NNMF
% 'snmf' -> Sparse NMF


% Algorithm = 'snmf';
% Sparsity = 0.3;

% Set number of synergies
MINK = 2;
MAXK = 2;

disp(['Estimated Time = ',num2str(((FinalSubject+1)*0.5)*(MAXK-MINK+1)), ' mins'])

% ------ Calibration phase START ---------- %

for Subject = StartSubject:FinalSubject
    
    [emg,force] = loadsubject(Subject);
    
    for k = MINK:MAXK
        disp(['Subject ',num2str(Subject), ' - k = ',num2str(k)]);
        
        Synergy = zeros(k*MAXDOF,electrodes);
        
        % DOF #1 (Index)    {0 - 50001}         3 repetitions
        % DOF #2 (Middle)   {100000 - 150000}   3 repetitions
        % DOF #3 (Ring)     {200000 - 250000}   3 repetitions
        % DOF #4 (little)   {295000 - 345000}   3 repetitions
        
        for DOF = 1:MAXDOF
            
            
            [start, finish] = loadindx(DOF,1);
            samples = finish - start;
            
            signal = emg(start:finish,1:electrodes);
            
            % Reescale emg signal before apply bayesian filtering
            option = 'Normalize';
            [Cal_data,Cal_force] = reescaling(signal',force',DOF,[],0,option);
            
            Cal_data = signal' * GAIN;
            
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
            
            if k==size(Z,1)
                W = ones(k,k);
                H = Z;
                fprintf('#SYN==#CHAN\n\t W = ones(%d,%d);\n\t H = Z;\n',k,k)
            else
                
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
        NinaResult(Subject).K(k).SYN = Synergy;
        NinaResult(Subject).K(k).CAL = [ControlDOF1;ControlDOF2;ControlDOF3];
        if MAXDOF == 4
            NinaResult(Subject).K(k).CAL = [ControlDOF1;ControlDOF2;ControlDOF3;ControlDOF4];
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
        figure;subplot(1,MAXDOF,1);hold on;plot(yd);plot(Csum); title('DOF #1')
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
        subplot(1,MAXDOF,2);hold on;plot(yd);plot(Csum);title('DOF #2')
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
        subplot(1,MAXDOF,3);hold on;plot(yd);plot(Csum);title('DOF #3')
        %axis([1 size(bayesSTD,2) maxRang])
        
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
            %axis([1 size(bayesSTD,2) maxRang])
        end
        disp('Done!')
        %close all
        
        % ------------- Regressor END ------------- %
        
        
        % -------------------- Calibration phase END ----------------------------- %
        
        
        %% --------------------- OFFLINE TEST ------------------------------- %%
        disp(['Offline Test start for Subject ', num2str(Subject)])
        
        % TASK #1 (Index)       {50000 - 100000} 3 reps
        % TASK #2 (Middle)      {165000 - 210000} 3 reps
        % TASK #3 (Ring)        {246000 - 296000} 3 reps
        % TASK #4 (Little)      {345000 - 395000} 3 reps
        % TASK #5 (Middle+Ring) {600000 - 650000} 3 reps
        % TASK #6 (Middle+Ring) {735000 - 780000} 3 reps
        
        
        % TASK = 1;
        
        for TASK = 1:6
            
            [Tstart, Tfinish] = loadindx(TASK,1);
            Tsamples = Tfinish - Tstart;
            Tsignal = emg(Tstart:Tfinish,1:electrodes);
            
            [Test_data,Test_force] = reescaling(Tsignal',force',DOF,[],0,option);
            
            Test_data = Tsignal' * GAIN; % Proportional scaling
            
            %figure
            %plot (Test_data)
            
            % -- Windowing parameters -- %
            % -- Perform bayessian filtering -- %
            
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
            
            NinaResult(Subject).K(k).TASK(TASK).RawCONT = Control;
            NinaResult(Subject).K(k).TASK(TASK).CONT = Sum;
            NinaResult(Subject).K(k).TASK(TASK).BAY = TbayesSTD;
            
            %             MVC_Force = max(max(For));
            %             Tol = 0.05; % Increasing the tolerance may decrease the euclidian distance
            %             Err = Tol*MVC_Force;
            
            disp(['Calculate errors for TASK #', num2str(TASK)]);
            for j = 1:4
                var1 = Sum(i,:);
                var2 = Tforce_res(:,i)';
                % var1(var1<0) = Err;
                % var2(var2<0) = Err;
                energy = (rms(var1).^2) * (2*size(var1,2)+1);
                NinaResult(Subject).K(k).ERR.D(TASK,j) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2));
                NinaResult(Subject).K(k).ERR.RMSE(TASK,j) = sqrt(sum((var1-var2).^2)/energy);
                NinaResult(Subject).K(k).ERR.RS(TASK,j) = Rsquare(var1, var2);
            end
            
        end
    end
disp(['ETA = ',num2str(((FinalSubject - Subject)*0.5)*(MAXK-MINK+1)), ' mins'])    
end
disp('Ouuuu yeah!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ///////////////////////////////////////////////

%% Plot Results for a specific Subject

% ----------------
Subject     =    37
TASK        =    1
k           =    2
% -----------------

figure('units', 'normalized', 'outerposition', [0 0.05 1 0.5])

for k = k:k
    [emg,force] = loadsubject(Subject);
    if Subject == 8
        Subject =9;
    end
    if TASK == 1
        Tstart = 50000;
        Tfinish = 100000;
    elseif TASK == 2
        Tstart = 150000;
        Tfinish = 200000;
    elseif TASK == 3
        Tstart = 250000;
        Tfinish = 300000;
    elseif TASK == 4
        Tstart = 345000;
        Tfinish = 395000;
    elseif TASK == 5 % DOF 1 and 4
        Tstart = 600000;
        Tfinish = 650000;
    elseif TASK == 6 % DOF 2 and 3
        Tstart = 730000;
        Tfinish = 780000;
    end
    
    % Reconstruct signal
    Sum1 = NinaResult(Subject).K(k).TASK(TASK).CONT(1,:);
    Sum2 = NinaResult(Subject).K(k).TASK(TASK).CONT(2,:);
    Sum3 = NinaResult(Subject).K(k).TASK(TASK).CONT(3,:);
    Sum = [Sum1;Sum2;Sum3];
    if MAXDOF == 4
        Sum4 = NinaResult(Subject).K(k).TASK(TASK).CONT(4,:);
        Sum = [Sum1;Sum2;Sum3;Sum4];
    end
    FOR = double(force);
    clear Tforce_res
    for i=1:MAXDOF
        Tforce_res(:,i) = ricampiona(FOR(Tstart:Tfinish,i),size(Sum,2),'linear');
    end
    clear S
    for i = 1:k*MAXDOF
        S(i,:) = NinaResult(Subject).K(k).SYN(i,:);
    end
    Range = [-5 max(max(max(Tforce_res),max(max(Sum))))+1];
    Range_Mid = [-5 max(max(max(Tforce_res),25))];
    lim = [0,11,1,20];
    samples = size(Sum,2);
    
    %DOF #1
    subplot(MAXDOF,6,k);hold on;plot(Tforce_res(:,1), 'LineWidth',2.5);plot(Sum1, 'LineWidth',2.5);
    axis([1 samples Range])
    title([ 'K = ',num2str(k)])
    if k == 1
        ylabel('Index')
    end
    subplot(MAXDOF,6,k+6);hold on;plot(Tforce_res(:,2),'LineWidth',2.5);plot(Sum2, 'LineWidth',2.5)
    axis([1 samples Range_Mid])
    if k+6 == 7
        ylabel('Middle')
    end
    subplot(MAXDOF,6,k+(6*2));hold on;plot(Tforce_res(:,3),'LineWidth',2.5);plot(Sum3, 'LineWidth',2.5)
    axis([1 samples Range])
    
    if k+(6*2) == 13
        ylabel('Ring')
    end
    if MAXDOF == 4
        subplot(MAXDOF,6,k+(6*3));hold on;plot(Tforce_res(:,4),'LineWidth',2.5);plot(Sum4, 'LineWidth',2.5)
        axis([1 samples Range])
        if k+(6*3) == 19
            ylabel('Little')
        end
    end
end


