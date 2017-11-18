clearvars -except Synergies Rawemg Comb1 Comb2 Comb3 Comb4 For1 For2 For3 For3 For4
% -----------------------------------------------
% Ideas for new adquisitions: 
% - Combination flexion/extension gestures
% - Simultaneous transition 0000-1000-1100-0100-0000
% - Complete ROM during calibration (flex to ext)
% -----------------------------------------------

% Design notch filter
% --------------------------------------------------------
Fs = 1200;
d  = fdesign.notch('N,F0,Q,Ap',30,462,100,1,Fs);
Hd = design(d);
%Filtemg = filter(Hd,emg);
% --------------------------------------------------------

Subject         =       1

% ------------- Calibration/Test data split ------------ %
% mode  = 3 ;  
% mode  = 2 ;   % All data
 mode  = 1 ;   % Select reps

REP = 4;

if REP == 1
Rep   = 1000:22000 ;
elseif REP == 2
Rep   = 22000:37000;
elseif REP == 3
Rep   = 37000:52000;
elseif REP == 4
Rep    = 1000:52000 ;
elseif REP == 12
Rep    = 1000:37000 ;
elseif REP == 23
Rep    = 22000:52000 ;
end

if Subject == 3 && REP == 3; Rep = 37000:49000; 
    elseif Subject == 3 && REP == 3; Rep = 1000:49000; 
        elseif Subject == 3 && REP == 23; Rep = 22000:49000;
             elseif Subject == 3 && REP == 4; Rep = 1000:49000;
end

Sdata  =    Rep    ; % Choose data for same trial calibration mode

Fdata  = (Sdata(1)/10):(Sdata(size(Sdata,2))/10);
% ------------------------------------------------------ %

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALIBRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% DOF-wise calibration: Estimate synergies of one finger at a time.
% ----------------------
Subject         =   Subject;
combination     =   'S' ; % Must be always 'S' during calibration
move            =   'f' ; % Flexion or extension
MAXDOF          =   4   ; % Maximun number of DOF's
k               =   2   ; % Number of synergies
trial           =   0   ; % Use different session during test ('0'
% ----------------------
trial_cal = trial;
% -- Perform calibration per each DOF -- % 


for DOF = 1:MAXDOF
    numberOffing    =   '1' ; % Must be always '1' during calibration
    force = [];
    % -- Load adquisition data -- %%
    finger = num2str(DOF+1);
    [emg,forces] = loadExoData(Subject,combination,numberOffing,finger,move,trial);
    force(DOF,:) = forces(str2double(finger)+1,:);   
    electrodes = size(emg,1)-1;
    %Filtemg = filter(Hd,emg);
    
    % -- Remove noise from the dataset 'mine_' -- %
    if Subject == 2 
        presignal = emg(2:end,1000:end); 
    else
        presignal = emg(2:end,100:end);
    end
    
    % -- Load calibration data -- %
    if mode == 1
        signal(:,:) = presignal(:,Sdata);
    elseif mode == 2
        signal = presignal;
    elseif mode == 3
        if DOF == 1
            signal = Comb1(2:end,:);
            force = For1;
        elseif DOF == 2
            signal = Comb2(2:end,:);
            force = For2;
        elseif DOF == 3
            signal = Comb3(2:end,:);
            force = For3;
        elseif DOF == 4
            signal = Comb4(2:end,:);
            force = For4;
        end
    end
    
    % -- Windowing parameters -- %
    samplingfreq = 1200 ;
    wTime        = 0.15 ; % Windows time in seconds
    over         = 50   ; % Windows overlap in percentage.

    % -- Reescale the emg signal before bayesian filter-- %
    option = 'Proportional'; 
    [Cal_data,Cal_force] = reescaling(signal,force,DOF,Fdata,mode,option);
    
    % -- Perform bayessian filtering -- %
    clear bayesSTD
    bins = 200;
    alpha = 10^-20;
    beta = 10^-500;
    MVC = 0.8;
    disp(['Bayes filter of DOF #', num2str(DOF)]);
    [bayesSTD force_res(DOF,:)] = windowing(Cal_data,Cal_force(DOF,:),wTime,samplingfreq,over,bins,alpha,beta,MVC);
  
    % ------------- NMF algorithm -------------
    disp(['Estimate Synergies of DOF #', num2str(DOF)]);
    % Start calibration
    Z = bayesSTD; % Each column of Z is a sample vector    
    % ---- Perform NNMF ---- %
    opt = statset('MaxIter',10); %'Display','final');
    [H0,W0] = nnmf(Z',k,'replicates',5,'options',opt,'algorithm','mult');
    opt = statset('Maxiter',100,'Display','final');
    [H,W] = nnmf(Z',k,'w0',H0,'h0',W0,'options',opt,'algorithm','als');
    H = H';
    W = W';
    % ---------------------- %
    
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
            Synergy(j+2*k,:) = W(:,j);
            ControlDOF3(j,:) = H(j,:);
        end       
    elseif DOF == 4
        ControlDOF4 = zeros(size(H));
        BayesDOF4 = bayesSTD;
        for j = 1:k
            Synergy(j+3*k,:) = W(:,j);
            ControlDOF4(j,:) = H(j,:);
        end
    end
    disp('Done!') 
    % ------------- NMF algorithm END -------------    
end
% -- Save the synergies for this session -- %
Synergies(Subject).REP(REP).K(k).SYN = Synergy;

disp('Calibration complete.')


% ----------- Scale factors  --------------- %

disp('Calculate regressor parameters...')
if move == 'e'
    maxRang = [-5 1];
else
    maxRang = [-1 8];
end
% DOF #1
yd = force_res(1,:)';
x_tr = [ControlDOF1(1,:)'];
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
if MAXDOF >= 2
    % DOF #2
    yd = force_res(2,:)';
    x_tr = [ControlDOF2(1,:)'];
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
end
if MAXDOF >= 3
    %DOF #3
    yd = force_res(3,:)';
    x_tr = [ControlDOF3(1,:)'];
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
end
% DOF #4
if MAXDOF == 4
    yd = force_res(4,:)';
    x_tr = [ControlDOF4(1,:)'];
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
%% ------------- Regressor END -------------

%%%%%%%%%%%%%%%%%%%%%%%%%%% CALIBRATION END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% OFFLINE TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FORCE_res force_res forces signal Test_force

% -------- Choose Fingers ------------- %
mode  = 1; % Select reps
Rep1  = 1000:22000 ;
Rep2  = 22000:37000;
Rep3  = 37000:52000;
All   = 1000:52000 ;

if Subject == 3
    Rep3 = 37000:49000; 
    All = 1000:49000 ;
end
Sdata_test = All; % One repetition for testing
Fdata_test = floor((Sdata_test(1)/10)):floor((Sdata_test(size(Sdata_test,2))/10));
if mode == 3; Fdata_test = 100; end

% ------------------------
trial           =    0    ;     % 0 for _I, 1 for _II, 2 for _III, etc
combination     =   'S'   ;     % Single or Multiple
DOF             =    1    ;     % Choose DOF (only if single):
fingers         =  [2 3 4]  ;     % Select simultaneous fingers
move            =   'f'   ;      % 'ef' special case
Classifier      =    1   ;
% ------------------------
trial_test = trial;
% ----------- Load Data ----------- %
if combination == 'S'
    numberOffing = '1';
    finger = num2str(DOF+1);
    [emg,forces] = loadExoData(Subject,combination,numberOffing,finger,move,trial);
elseif combination == 'M'
    numberOffing = char(string((size(fingers,2))));
    finger = num2str(fingers(1));
    for j = 2:str2double(numberOffing)
        finger = strcat(finger,num2str(fingers(j)));
    end
    [emg,forces] = loadExoData(Subject,combination,numberOffing,finger,move,trial);
end
% --------- Load Data END --------- %


% --- Preprocess data before filtering --- %
electrodes = size(emg,1)-1; % May be different for each subject
samples = size(emg,2);
emg = filter(Hd,emg); % Notch filter at 460 Hz
presignal = emg(2:end,100:end);
if Subject == 2 ;presignal = emg(2:end,1000:end); end
if mode == 1
    signal(:,:) = presignal(:,Sdata_test);
else
    signal = presignal;
end

% -- Reescale the signal -- %
[Test_data,Test_force] = reescaling(signal,forces,DOF,Fdata_test,mode,option);

% -- Bayesian filter -- %
clear bayesSTD
disp('Bayes filter...');
[bayesSTD force_res(DOF,:)] = windowing(Test_data,Test_force(DOF,:),wTime,samplingfreq,over,bins,alpha,beta,MVC);
% -- Bayesian filter END -- %

% ------ NMF -------- %
disp('Estimating Control signals...');
Z = bayesSTD; % Each column of Z is a sample vector
SynInv = pinv(Synergy)';
%Control = NMFwindowing(Z,SynInv,wTime,samplingfreq,over);
Control = SynInv * Z; % C[mxp] = A[2*MAXDOF x channels] * B[channels x samples] Use Moore-Penrose algorithm.
disp('Done!');
% ------ NMF END -------- %


% ------ Resampling forces ------ %
% Forces(1)=time; (2)=thumb; (3)=Index; (4)=middle; (5)=ring; (6)=little
if mode == 1
    FOR = double(forces(3:MAXDOF+2,Fdata_test));
elseif mode == 2
    FOR = double(forces(3:MAXDOF+2,:));
end
for j = 1:size(FOR,1)
    FORCE_res(j,:) = ricampiona(FOR(j,:),size(Z,2),'linear');
end
% -------------------------------- %

% /////// Classifier Test ///////
if Classifier == 1
    data = 0; % data = 0 for exoskeleton data
    for s = 1:MAXDOF
        [REstBayes,Rbayes,Terror,Fing] = classifier(bayesSTD,Control,Synergy,s,0,k,data);
        Terror
        if Terror > 15
           disp(['DOF : ',num2str(s),' Inactive'])
            Control(Fing,:) = Control(Fing,:)*0.1;
        end
        EstBay{s} = REstBayes;
    end
end
% ////////////////////////////////////////////

% ------- Reconstruct Signal ------- %
clear C

for i = 1:k*MAXDOF
    C(i,:) = Control(i,:);
end
ca1 = C(1,:)*beta1_1(2);
ca2 = C(k+1,:)*beta2_1(2);
ca3 = C(2*k+1,:)*beta3_1(2);
ca4 = C(3*k+1,:)*beta4_1(2);
for i = 2:k
    ca1 = ca1 + (C(i,:) * beta1_1(i+1));
end
for i = 2:k
    ca2 = ca2 + (C(i+k,:) * beta2_1(i+1));
end
for i = 2:k
    ca3 = ca3 + (C(i+k*2,:) * beta3_1(i+1));
end
for i = 2:k
    ca4 = ca4 + (C(i+k*3,:) * beta4_1(i+1));
end
Sum1 = beta1_1(1) + ca1;
Sum2 = beta2_1(1) + ca2;
Sum3 = beta3_1(1) + ca3;
Sum4 = beta4_1(1) + ca4;
% -- Nonnegative control signal -- %
% disp('Absolute values active!!');
% Sum1 = abs(Sum1);
% Sum2 = abs(Sum2);
% Sum3 = abs(Sum3);
% Sum4 = abs(Sum4);

% -------------------------------- %
% ------- Plot results ------- %
if move == 'e'
    maxRang = [-5 1];
else
    maxRang = [-1 8];
end

figure;
subplot(MAXDOF,1,1);hold on;plot(FORCE_res(1,:));plot(Sum1);
axis([1 size(bayesSTD,2) maxRang])
title(['\fontsize{12}Subject ',num2str(Subject),' - DOF(s) ',finger, ...
    ' active - Calibration data Rep: ',num2str(REP),' & Trial: ',num2str(trial_cal)])
ylabel('Index')
subplot(MAXDOF,1,2);hold on;plot(FORCE_res(2,:));plot(Sum2)
axis([1 size(bayesSTD,2) maxRang])
ylabel('Middle')
subplot(MAXDOF,1,3);hold on;plot(FORCE_res(3,:));plot(Sum3)
axis([1 size(bayesSTD,2) maxRang])
ylabel('Ring')
subplot(MAXDOF,1,4);hold on;plot(FORCE_res(4,:));plot(Sum4)
axis([1 size(bayesSTD,2) maxRang])
ylabel('Little')
xlabel(['Test data Trial: ',num2str(trial_test)])
%%%%%%%%%%%%%%%%%%%%%%%%%%% OFFLINE TEST END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%