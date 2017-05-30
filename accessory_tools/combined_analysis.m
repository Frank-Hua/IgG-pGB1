%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

This is for combining individual output files together and analyzing them.

Check and adjust parameters that are marked with "frank".
%}

function combined_analysis

%initial input
path=command_input('input directory','C:\\Users\\frank\\Documents\\MATLAB','s');
cd([path '\molecules']);

A=dir;
[nf,dum]=size(A);

%{
countsT generates intensity histograms of each molecule and makes sure
    each molecule contributes equally.
traces gather all traces in one variable.

%}
% molecules=[];
% intensities=[];
countsT=[];
traces=[];

xbins=-50:25:2500;
for i=1:nf,
    if A(i).isdir == 0
        s=A(i).name;
        if strcmp(s(end-2:end), 'txt')
            disp(s);
%             cd([path '\molecules']);
%             temp=dlmread(s,'\t');
%             molecules=[molecules,temp];
            
            cd([path '\intensities']);
            s1=['intensity' s(9:end)];
            intensity=dlmread(s1,'\t');
            [counts,centers] = hist(intensity,xbins);
%             f = fit(centers.',counts.','gauss1');
%             intensities=[intensities,f.b1];
            counts=counts/sum(counts);
            countsT=[countsT,counts'];
            
            cd([path '\traces']);
            s2=['trace' s(9:end)];
            temp2=dlmread(s2,'\t');
            traces=[traces,temp2];
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
fn='intensity_histogram.txt';
save(fn,'countsT','-ascii');

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
tr2(1:100)=0;
tr_label=bwlabel(tr2);
props=regionprops(tr_label,'Area');
temp_off=[props.Area];

end