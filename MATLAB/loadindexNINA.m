function [s,f] = loadindexNINA(D,mod)

% Load indexes for NinaPro Database V2.
% mod = 0 -> Calibration data. 3 first reps. Only single DOF's. 4 Tasks
% mod = 1 -> Testing data. 3 last reps. Single and simultaneous fingers. 6 Tasks


if mod == 0
    
    if D == 1
        s = 1;
        f = 50001;
    elseif D == 2
        s = 100000;
        f = 150000;
    elseif D == 3
        s = 200000;
        f = 250000;
    elseif D == 4
        s = 295000;
        f = 345000;
    end
    
elseif mod == 1
    
    if D == 1
        s = 50000;
        f = 100000;
    elseif D == 2
        s = 150000;
        f = 200000;
    elseif D == 3
        s = 250000;
        f = 300000;
    elseif D == 4
        s = 345000;
        f = 395000;
        % DOF 1 and 4
    elseif D == 5
        s = 600000;
        f = 650000;
        % DOF 2 and 3
    elseif D == 6
        s = 730000;
        f = 780000;
    end
end
end