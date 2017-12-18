%% Choose Subject and Test DOF
% ----------------
Subject     =    1
TASK        =    1
k           =    2
Classifier  =    0
% -----------------
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% /////////// Classifier test ////////////////////
NinaResult = Subjects40_4DOFS;
Bastards = [];
TotalError = zeros(MAXK,MAXDOF);
% for k = 1:6
%     disp([' Test of k = ',num2str(k)])
% for TASK = 1:MAXDOF
    Match = 0;
    fail = 0;
    for Subject = 1:20
        
[emg,force] = loadsubject(Subject);
    if Subject == 8
        Subject = 9;
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
if MAXDOF == 4
    Sum4 = NinaResult(Subject).K(k).TASK(TASK).CONT(4,:);
end
FOR = double(force);
clear Tforce_res
for i=1:MAXDOF
    Tforce_res(:,i) = ricampiona(FOR(Tstart:Tfinish,i),size(Sum1,2),'linear');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /////////////// Classifier ///////////////////
Syn = NinaResult(Subject).K(k).SYN;
Bayess = NinaResult(Subject).K(k).TASK(TASK).BAY;
RawCont = NinaResult(Subject).K(k).TASK(TASK).RawCONT;
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
    if Active == TASK
        Match = Match+1;
    else
       %disp(['Fail!! of Subject ',num2str(Subject),' for DOF: ',num2str(TASK),', K = ', num2str(k)])
       fail = [fail Subject];
    end

    end
    ClasErr(TASK) = Match;
    
    if TASK == 1
        Bastards1 = fail;
    elseif TASK == 2
        Bastards2 = fail;
    elseif TASK == 3
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
%% Plot errors
Subject = 37
close all
MINK = 1;
MAXK = 6;
clear TD1 TD2 TD3 TD4 TD5 TD6
for i = MINK:MAXK
    TD1(i,:) = NinaResult(Subject).K(i).ERR.RMSE(1,:);
    TD2(i,:) = NinaResult(Subject).K(i).ERR.RMSE(2,:);
    TD3(i,:) = NinaResult(Subject).K(i).ERR.RMSE(3,:);
    TD4(i,:) = NinaResult(Subject).K(i).ERR.RMSE(6,:);
    TD = max([TD1 TD2 TD3 TD4]);
    if MAXDOF == 4
    TD4(i,:) = NinaResult(Subject).K(i).ERR.RMSE(4,:);
    TD5(i,:) = NinaResult(Subject).K(i).ERR.RMSE(5,:);
    TD6(i,:) = NinaResult(Subject).K(i).ERR.RMSE(6,:);
    TD = max([TD1 TD2 TD3 TD4 TD5 TD6]);
    end
end

subplot(2,3,1)
bar(TD1./size(NinaResult(Subject).K(1).CAL,2))
title('Error for DOF #1 active')
xlabel('k')
subplot(2,3,2)
bar(TD2./size(NinaResult(Subject).K(1).CAL,2))
title('Error for DOF #2 active')
subplot(2,3,3)
bar(TD3./size(NinaResult(Subject).K(1).CAL,2))
title('Error for DOF #3 active')
subplot(2,3,4)
bar(TD4./size(NinaResult(Subject).K(1).CAL,2))
title('Error for DOFs #4 active')
if MAXDOF == 4
subplot(2,3,5)
bar(TD5./size(NinaResult(Subject).K(1).CAL,2))
title('Error for DOFs #1 and #4 active')
subplot(2,3,6)
bar(TD6./size(NinaResult(Subject).K(1).CAL,2))
title('Error for DOFs #2 and #3 active')
end
%bar(NinaResult(Subject).SYN(1,:));hold all;bar(NinaResult(2).SYN(1,:));
%% -- Plot error of active DOFs (1 to 4) -- %%
clear E1 E2 E3 E4 E5 E6 
for m = 1:FinalSubject
    for n = 1:MAXDOF
    E1(m,n) = NinaResult(m).K(1).ERR.RS(n,n);
    end
end

for m = 1:FinalSubject
    for n = 1:MAXDOF
    E2(m,n) = NinaResult(m).K(2).ERR.RS(n,n);
    end
end
for m = 1:FinalSubject
    for n = 1:MAXDOF
    E3(m,n) = NinaResult(m).K(3).ERR.RS(n,n);
    end
end
for m = 1:FinalSubject
    for n = 1:MAXDOF
    E4(m,n) = NinaResult(m).K(4).ERR.RS(n,n);
    end
end
for m = 1:FinalSubject
    for n = 1:MAXDOF
    E5(m,n) = NinaResult(m).K(5).ERR.RS(n,n);
    end
end
for m = 1:FinalSubject
    for n = 1:MAXDOF
    E6(m,n) = NinaResult(m).K(6).ERR.RS(n,n);
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