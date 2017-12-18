function R_Square = Rsquare(estimates, targets)

% Function to compute the R² index proposed by d'Avella
%
% This value will represent one value for all DOF: for the whole
% estimation
%
% 'estimates' is expected to be a matrix(#DOF, #Samples), it contains the
% values estimated by the model
% 'targets' is expected to be a matrix(#DOF, #Samples), it contains the
% measured target values

%  It provides the percentage of total variation in the force functions (targets)captured by the
% estimates. The larger the value, the better the estimation performance. 
% The index assumes its upper bound of unity in the case
% of perfect estimate (zero mse). Its lower bound is -inf (evil!!!!)

if ((size(estimates,1)~=size(targets,1))||(size(estimates,2)~=size(targets,2)))
    error('Estimates and targets have not the same size')
end
nDoF = size(estimates, 1);
nSamples = size(estimates, 2);
avgTargets = mean(targets, 2); %column vector [nDof x1]
avgTargetsMatr = avgTargets .*ones(1,size(targets,2)); %[nDof x #Samples]

for DoF_ix = 1:nDoF
clear numerator denominator
    numerator   = (sum( (estimates(DoF_ix,:)-    targets(DoF_ix,:)    ).^2));%SSres
    denominator = (sum( (targets(DoF_ix,:)  -   avgTargetsMatr(DoF_ix,:)).^2));%SStot
    R_Square(DoF_ix) = 1 - (numerator ./ denominator);
end


% nDoF = size(estimates, 1);
% nSamples = size(estimates, 2);
% numerator = 0;
% denominator =0;
% for DoF_ix = 1:nDoF
%     numerator = numerator+ sum((estimates(DoF_ix,:) - targets(DoF_ix,:)).^2);%,2);
%     denominator = denominator + sum((targets(DoF_ix,:) - avgTargets(DoF_ix)*ones(1,nSamples)).^2);%,2);
% end   
% R_Square = 1 - (numerator ./ denominator);


end



%%%% Funzione Che usava Francesco per stimare la r di ricostruzione delle
%%%% sinergie
%%%  Ricorda di sostituire   
% % %     dataMat = targets;
% % %     W*H = estimates;

% SSE=sum(sum(       (dataMat'-W*H).^2          ));
% for i=1:LExos.nEMGsignals
% 	SST(i,:)=detrend(dataMat(i,:),'constant').^2;
% end
% SST=sum(sum(SST));
% rSquared=1- SSE/SST ;


%%%% fatte le sostituzioni
% SSE=sum(sum(       (targets-estimates).^2          ));
% for i=1:LExos.nEMGsignals
% 	SST(i,:)=detrend(targets(i,:),'constant').^2;
% end
% SST=sum(sum(SST));
% rSquared=1- SSE/SST ;