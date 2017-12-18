function [c d f]= loadTask(T)

if T == 1
    c    =   'S'   ;     % Single or Multiple
    d    =    1    ;     % Choose DOF (only if single):
    f    =    0    ;
elseif T == 2
    c    =   'S'   ;     % Single or Multiple
    d    =    2    ;     % Choose DOF (only if single):
    f    =    0    ;
elseif T == 3
    c    =   'S'   ;     % Single or Multiple
    d    =    3    ;     % Choose DOF (only if single):
    f    =    0    ;
elseif T == 4
    c    =   'S'   ;     % Single or Multiple
    d    =    4    ;     % Choose DOF (only if single):
    f    =    0    ;
elseif T == 5
    d    =    1    ;
    c    =   'M'   ;     % Single or Multiple     % Choose DOF (only if single):
    f    =  [2 3]  ;     % Select simultaneous fingers
elseif T == 6
    d    =    1    ;
    c    =   'M'   ;     % Single or Multiple     % Choose DOF (only if single):
    f    =  [3 4]  ;     % Select simultaneous fingers
elseif T == 7
    d    =    1    ;
    c    =   'M'   ;     % Single or Multiple     % Choose DOF (only if single):
    f    =  [4 5]  ;
elseif T == 8
    d    =    1    ;
    c    =   'M'   ;     % Single or Multiple     % Choose DOF (only if single):
    f    =  [2 3 4]  ;
elseif T == 9
    d    =    1    ;
    c    =   'M'   ;     % Single or Multiple     % Choose DOF (only if single):
    f    =  [3 4 5]  ;
elseif T == 10
    d    =    1    ;
    c    =   'M'   ;     % Single or Multiple     % Choose DOF (only if single):
    f    =  [2 3 4 5]  ;
end
