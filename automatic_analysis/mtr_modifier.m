%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function clicktime=mtr_modifier(index,len,clicktime)

[X,~,button] = ginput;
if isempty(X) || mod(length(X),2) ~= 0
    return;
end

X=round(X);

index2 = find(button==1);
index3 = find(button==3);

X2=X(index2);
temp = find(X2<1);
X2(temp) = 1;
temp = find(X2>len);
X2(temp) = len;

X3=X(index3);
temp = find(X3<1);
X3(temp) = 1;
temp = find(X3>length(index));
X3(temp) = length(index);

clicktime = [clicktime' X2']';
clicktime = [clicktime' index(X3)']';

return;
end