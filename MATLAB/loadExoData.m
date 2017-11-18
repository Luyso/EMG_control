% Subject = 1 to ...
% combination: 'S' single or 'M' multiple
% numberOffing: 1 , 2
% finger: '1' thumb, '2' index, '3' middle, '4' ring, '5' little
% move: 'f' flexion or 'e' extension.

function [emg,force] = loadExoData(Subject,combination,numberOffing,finger,move,test)
if Subject == 1
    name = 'Luis_';
elseif Subject == 2
    name = 'Mine_';
elseif Subject == 3
    name = 'Louis_';
elseif Subject == 4
    name = 'Luis_Bipolar_';
end

if test == 1
    rep = '_II.mat';
elseif test == 2
    rep = '_III.mat';
else
    rep = '_I.mat';
end

d = strcat(name,combination,'_',numberOffing,'_',finger,'_',move,rep);

load(d,'RawEMG','EMG','Forces')

emg = EMG;
force = Forces;
end