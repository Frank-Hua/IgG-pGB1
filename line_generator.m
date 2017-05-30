%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function [line_x,line_y]=line_generator(index)

d_tr=conv(index,[1,-1],'same');
d_tr(end)=1;
index2=find(d_tr>1);

line_x=[index2 index2]';
line_y=[zeros(length(index2),1) ones(length(index2),1)]';

end