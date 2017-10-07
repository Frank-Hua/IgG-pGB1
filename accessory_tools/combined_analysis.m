%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

This is for combining individual trace files together and analyzing them.

Adjustable parameters are marked with "frank".
%}

function combined_analysis

%% initial input
path=command_input('input directory:','C:\\Users\\frank\\Documents\\MATLAB','s');
cd([path '\traces']);

A=dir;
[nf,~]=size(A);

%{
traces gather all traces together.
%}

traces=[];

for i=1:nf
    if A(i).isdir == 0
        s=A(i).name;
        if strcmp(s(end-2:end), 'txt')
            disp(s);
            temp=dlmread(s,'\t');
            traces=[traces,temp];
        end
    end
end

len=size(traces,1);
n=size(traces,2);

on_durations=[];
off_durations=[];
for i=1:n
    [temp_on,temp_off]=durations(traces(:,i)',len);
    on_durations=[on_durations,temp_on];
    off_durations=[off_durations,temp_off];
end

cd(path);

fn='on_durations.txt';
fid=fopen(fn,'w');
fprintf(fid, '%f\n', on_durations);
fclose(fid);

fn='off_durations.txt';
fid=fopen(fn,'w');
fprintf(fid, '%f\n', off_durations);
fclose(fid);

close('all');
end







function [temp_on,temp_off]=durations(tr,len)

tr_label=bwlabel(tr);
props=regionprops(tr_label,'Area');
temp_on=[props.Area];

tr2=ones(1,len)-tr;
tr_label=bwlabel(tr2);
props=regionprops(tr_label,'Area');
temp_off=[props.Area];

end