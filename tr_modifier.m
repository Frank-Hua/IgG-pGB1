%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Adjustable parameters are marked with "frank".
%}

function [clicktime2,clicktime3]=tr_modifier(index,clicktime2,clicktime3)

input('enter-to continue ','s');

[X,~,button] = ginput;
if isempty(X) || mod(length(X),2) ~= 0
    return;
end

X=round(X);

temp = find(X<1);
X(temp) = 1;
temp = find(X>length(index));
X(temp) = length(index);

index2 = find(button==1);
index3 = find(button==3);
clicktime2 = [clicktime2' index(X(index2))']';
clicktime3 = [clicktime3' index(X(index3))']';

return;
end