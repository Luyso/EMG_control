function [BAYES, FORCE_res] = bayes(data,FORCE,wTime,fs,over)
%% clear all; clc; close all;
% TT = DATA_ALL(1,cursor_info(2).DataIndex:cursor_info(1).DataIndex);
% EMG = DATA_ALL([3:4],cursor_info(2).DataIndex:cursor_info(1).DataIndex);
% FORCE = DATA_ALL(5,cursor_info(2).DataIndex:cursor_info(1).DataIndex);
% REFERENCE= DATA_ALL(end,cursor_info(2).DataIndex:cursor_info(1).DataIndex);
%load('EMG_PROVA_withFORCE.mat')
EMG = data'; %[1 x Samples]
%EMG = signal(1:50000,1)';
%fea = 'RMS';
%EMG = Z(1,:);
%FORCE = force(1:50000,1)';

%EMG(2,:) = [];
% figure
% plot_array([EMG;FORCE],{'EMG','FORCE'})
% plot_array([EMG;FORCE])
%FS = round(1/mean(diff(TT)));
FS = fs;
% estrazione features
WindowSize_ms = wTime; %0.100; %seconds
WindowSize = floor((WindowSize_ms * FS) / 5 );
OVERLAP = over; %percentage
WindowIncrement = WindowSize - round(WindowSize*OVERLAP/100);
LengthData = size(EMG,2);
LENGTH_FEAT = floor((LengthData - WindowSize )/WindowIncrement + 1);

%
CurrIx = 1;
BAYf = [];
% -------------  Parameters for Bayesian filter  -------------
    param.sf = fs;  % set sampling frequency
    param.bins = 500;         % output i.e. parameter quantization levels (80)
    param.alpha = 10^-9;     % sets diffusion rate (-40)
    param.beta = 10^-900;      % sets probability of sudden jumps (-100)
    param.model = 'Laplace';  % choose between 'Laplace' and 'Gauss'
    param.pointmax = false;   % false: use expectation value as point estimation,
    % true: use maximum of posterior as point estimation
    param.sigmaMVC = 1;     % MVC amplitude value (80% for NinaPro DB2)
    % -------------  Parameters for Bayesian filter END  -------------
    pri = ones(param.bins,1)/param.bins;         % define uniform prior
    x = EMG/param.sigmaMVC;   % rescale data with respect to sigma MVC. This
    % helps avoiding numerical problems in case of the raw data being
    % especially if not single values are processed with the function
    % BayesFilter but instead multiple EMG measurements at a time. Refer also
    % to the comments in BayesFilter.m.
    %
    % plot(x);
    %
    sig = linspace(param.bins^-1,1,param.bins)';  % sigma (amplitude) axis can
    % start at 0 or param.bins^-1. This is a matter of taste. Note however that
    % sigma=0 will always have zero probability by definition (0 is an absorbing
    % boundary) due to the divergence at 0 when computing the likelihood. Refer
    % also to the comments in BayesFilter.m.
    %
    % perform the filtering
while (CurrIx+WindowSize-1) <= LengthData

    CurrWindow = CurrIx:WindowSize+CurrIx-1;
    temp_data = x(:,CurrWindow);
    temp_bay = zeros(size(temp_data));
%       for j = 1:size(temp_data,2)
%         for i = 1:size(temp_data,1)
              [temp_bay, pri] = BayesFilter(temp_data, pri, sig, param);
%         end
%       end
    BAYf     = [BAYf, temp_bay];
    CurrIx = CurrIx+WindowIncrement;
    
end
%
%
% FEAT = [RMSf;MAVf;ZCf;SSCf;WLf;VAR;VARlog];
% figure
%plot_array(FEAT,{'RMS','MAV','ZC','SSC','WL','Var','logVar'});
%plot_array(BAYf);
% %
%display(['Feature choosen: ',fea])
BAYES = BAYf;
%
% %%  RESAMPLING THE FORCE IN ORDER TO MATCH THE FEATURES
% % FORCE_res2 = ricampiona(FORCE,length(FEAT));
 FOR = double(FORCE);
 %FORCE_res = resample(FOR,LENGTH_FEAT,LengthData);
 FORCE_res = ricampiona(FOR,LENGTH_FEAT,'linear');
% figure
% subplot(4,1,1)
% plot(RMSf);
% subplot(4,1,2)
% plot(FORCE_res)
% subplot(4,1,3)
% plot(EMG_ACT)
% subplot(4,1,4)
% plot(FORCE)
% 
% legend('using resample','using "ricampiona"')
% 
% 
end
% %%
% LRmyOBJ = LR_WithReg()
% LRmyOBJ.USEBIAS = 0;
% LRmyOBJ.useRegularization = 0;
% LRmyOBJ.train(FEAT,FORCE_res)
% %%
% Y_hat = LRmyOBJ.apply(FEAT)
% % Y_hat = W'*FEAT;
% LRmyOBJ.showModel
% figure
% % imagesc(W)
% figure
% plot(FORCE_res);hold on;plot(Y_hat)
% legend('Measured Force - downsampled','estimated force')
% 
