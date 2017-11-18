function [emg,force] = loadsubject(s)

Subject = s;

if Subject == 1
    load('S1_E3_A1.mat','emg','force')
elseif Subject == 2
    load('S2_E3_A1.mat','emg','force')
elseif Subject == 3
    load('S3_E3_A1.mat','emg','force')
elseif Subject == 4
    load('S4_E3_A1.mat','emg','force')
elseif Subject == 5
    load('S5_E3_A1.mat','emg','force')
elseif Subject == 6 % EMG not good
    load('S6_E3_A1.mat','emg','force')
elseif Subject == 7 % EMG not good. Forces not good
    load('S7_E3_A1.mat','emg','force')
elseif Subject == 8
    load('S8_E3_A1.mat','emg','force')
elseif Subject == 9 % Force really bad
    load('S9_E3_A1.mat','emg','force')
elseif Subject == 10
    load('S10_E3_A1.mat','emg','force')
elseif Subject == 11 % Force really bad
    load('S11_E3_A1.mat','emg','force')
elseif Subject == 12 % Seems ok
    load('S12_E3_A1.mat','emg','force')
elseif Subject == 13
    load('S13_E3_A1.mat','emg','force')
elseif Subject == 14 % Looks nice
    load('S14_E3_A1.mat','emg','force')
elseif Subject == 15 % Not bad
    load('S15_E3_A1.mat','emg','force')
elseif Subject == 16
    load('S16_E3_A1.mat','emg','force')
elseif Subject == 17
    load('S17_E3_A1.mat','emg','force')
elseif Subject == 18
    load('S18_E3_A1.mat','emg','force')
elseif Subject == 19
    load('S19_E3_A1.mat','emg','force')
elseif Subject == 20
    load('S20_E3_A1.mat','emg','force')
elseif Subject == 21
    load('S21_E3_A1.mat','emg','force')
elseif Subject == 22
    load('S22_E3_A1.mat','emg','force')
elseif Subject == 23
    load('S23_E3_A1.mat','emg','force')
elseif Subject == 24
    load('S24_E3_A1.mat','emg','force')
elseif Subject == 25
    load('S25_E3_A1.mat','emg','force')
elseif Subject == 26
    load('S26_E3_A1.mat','emg','force')
elseif Subject == 27
    load('S27_E3_A1.mat','emg','force')
elseif Subject == 28
    load('S28_E3_A1.mat','emg','force')
elseif Subject == 29
    load('S29_E3_A1.mat','emg','force')
elseif Subject == 30
    load('S30_E3_A1.mat','emg','force')
elseif Subject == 31
    load('S31_E3_A1.mat','emg','force')
elseif Subject == 32
    load('S32_E3_A1.mat','emg','force')
elseif Subject == 33
    load('S33_E3_A1.mat','emg','force')
elseif Subject == 34
    load('S34_E3_A1.mat','emg','force')
elseif Subject == 35
    load('S35_E3_A1.mat','emg','force')
elseif Subject == 36
    load('S36_E3_A1.mat','emg','force')
elseif Subject == 37
    load('S37_E3_A1.mat','emg','force')
elseif Subject == 38
    load('S38_E3_A1.mat','emg','force')
elseif Subject == 39
    load('S39_E3_A1.mat','emg','force')
elseif Subject == 40
    load('S40_E3_A1.mat','emg','force')
end

emg = emg;
force = force;

end
