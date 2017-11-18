function [REstBayes,Rbayes, Terror,Fing] = classifier(bayess,Control,Synergy,s,plots,k,data)
        if s == 1
            Fing = 1:k;
        elseif s==2
            Fing = k+1:2*k;
        elseif s == 3
            Fing = 2*k+1:3*k;
        elseif s == 4
            Fing = 3*k+1:4*k;
        end
        
        EstBayes = Synergy(Fing,:)' * (Control(Fing,:));
        % -- Rescaling -- %%
        Rbayes = bayess / max(max(bayess));
        REstBayes = EstBayes + abs(min(min(EstBayes)));
        REstBayes = REstBayes / max(max(REstBayes));
        % -- Normalize -- %%
        m = mean(mean(bayess));
        d = std(std(bayess));
        Rbayes = (bayess - m) / d;
        m = mean(mean(EstBayes));
        d = std(std(EstBayes));
        REstBayes = (EstBayes - m) / d;
        
        if data == 1 % If used with ninapro mode == 1 for Exohand mode == 0;
            Rbayes = bayess;
            REstBayes = EstBayes;
        end
        
        % --------------- %%
        %dist = dtw(Rbayes,REstBayes');
        Terror = 0;
        Teuc = 0;
        Error = [];
        for m = 1:size(Rbayes,1)
            var1 = Rbayes(m,:);
            var2 = REstBayes(m,:);
%             var1(var1<0) = 0.1;
%             var2(var2<0) = 0.1;
%             Euc = (var1-var2)' * (var1-var2);
%             Teuc = Teuc + Euc;
            Error(m) = dot(var1-var2, var1-var2)/sqrt(dot(var1,var1)*dot(var2,var2));
        end
%        disp(['DOF ',num2str(s),':'])
        Terror = sum(Error); %Terror + Error;
        Merror = mean(Error)
        Error;
        if plots == 1
        figure
        ax = subplot(1,1,1);
        plot(REstBayes)
        title(['Estimated bayesian for DOF: ',num2str(s),'. Euclidian distance = ',num2str(Terror)])
        ax.XTick = [];
        figure;
        ax = subplot(1,1,1);
        plot(Rbayes','linewidth',1);
        ax.XTick = [];
        title('Bayesian signal')
        end
%         if Terror > 15
%             Control(Fing,:) = Control(Fing,:)*0.1;
%         end
    
    

end
 