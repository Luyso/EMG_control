% Questa funzione prende un vettore (input) e lo ricampiona sulla base del
% numero di campioni richiesti (n)
% output=ricampiona(input,n)
% method = interpolazione ('linear','cubic','spline')

function output=ricampiona(input,n,method)


c=length(input);
asc1=1:c;
asc2=0:c/n:c;

if nargin<3
    disp('Spline Interpolation')
    output=spline(asc1,input,asc2);
else
    output=interp1(asc1,input,asc2,method);
end



output(1)=[];
