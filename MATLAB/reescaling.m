function [Cal_dat,Cal_forc] = reescaling(sig,forc,DOF,Fdata,mod,opt)

Nsignal = zeros(size(sig));
Rsignal = zeros(size(sig));
RNsignal = zeros(size(sig));
NRsignal = zeros(size(sig));
NPsignal = zeros(size(sig));
Psignal = zeros(size(sig));
Msignal = zeros(size(sig));
for i = 1:size(sig,1)
        % Normalization
        m = mean(sig(i,:));
        d = std(sig(i,:));
        Nsignal(i,:) = (sig(i,:) - m) / d;
        % Rescaling
        Rsignal(i,:) = sig(i,:) - min(sig(i,:));
        Rsignal(i,:) = (Rsignal(i,:)/range(Rsignal(i,:)))*(2);
        Rsignal(i,:) = (Rsignal(i,:) -1);
        % Rescaling + Normalization
        m = mean(Rsignal(i,:));
        d = std(Rsignal(i,:));
        RNsignal(i,:) = (Rsignal(i,:) - m) / d;
        % Normalization + Rescaling
        NRsignal(i,:) = Nsignal(i,:) - min(Nsignal(i,:));
        NRsignal(i,:) = (NRsignal(i,:)/range(NRsignal(i,:)))*(1 - (-1));
        NRsignal(i,:) = (NRsignal(i,:) - 1);
        % Normalization + proportional rescaling (/max)
        NPsignal(i,:) = Nsignal(i,:) / (max(max(Nsignal(i,:))));
        % Proportional
        Psignal(i,:) = sig(i,:) / 400; 
        % MAX
        Msignal(i,:) = sig(i,:) / (max(max(sig(i,:))));
        % Choose the pre-processed signal
        if     strncmp(opt,'Reescale',size(opt,2))
            Cal_dat(i,:) = Rsignal(i,:);
        elseif strncmp(opt,'Normalize',size(opt,2))
            Cal_dat(i,:) = NPsignal(i,:);
        elseif strncmp(opt,'Proportional',size(opt,2))
            Cal_dat(i,:) = Psignal(i,:);
        elseif strncmp(opt,'Max',size(opt,2))
            Cal_dat(i,:) = Msignal(i,:);
        else
            disp('Reescaling option unknown')
        end
        % When mode = 1 use only choosen repetition for calibration
        if mod == 1
            Cal_forc(DOF,:) = forc(DOF,Fdata);
        else
            Cal_forc(DOF,:) = forc(DOF,:);
        end
        
end
end
