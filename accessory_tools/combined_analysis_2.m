%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

This is for combining individual output files and analyzing them together.
This one has a threshold, as compared to just combined_analysis.m

Check and adjust parameters that are marked with "frank".
%}

function combined_analysis_2

%initial input
path=command_input('input directory','C:\\Users\\frank\\Documents\\MATLAB','s');
cd([path '\molecules']);

A=dir;
[nf,dum]=size(A);

%{
molecules gather all molecules' info together in one file.
intensities store intensity peak positions of each molecule.
countsT generates intensity histograms for each molecule and makes sure
    each molecule contributes equally.
traces gather all traces in one file.

%}
molecules=[];
intensities=[];
countsT=[];
traces=[];

xbins=-50:25:2500;
n=0;
%thres=0.0017*15100;
thres=0.0008*15100;
for i=1:nf,
    if A(i).isdir == 0
        s=A(i).name;
        if strcmp(s(end-2:end), 'txt')
            disp(s);
            cd([path '\molecules']);
            temp=dlmread(s,'\t');
            if temp(3) <= thres;
                molecules=[molecules,temp];

                cd([path '\intensities']);
                s1=['intensity' s(9:end)];
                intensity=dlmread(s1,'\t');
                [counts,centers] = hist(intensity,xbins);
                f = fit(centers.',counts.','gauss1');
                intensities=[intensities,f.b1];
                counts=counts/sum(counts);
                countsT=[countsT,counts'];

                cd([path '\traces']);
                s2=['trace' s(9:end)];
                temp2=dlmread(s2,'\t');
                traces=[traces,temp2];
                n=n+1;
            end
        end
    end
end

disp(mean(molecules(1,:)));
disp(std(molecules(1,:)));
disp(size(molecules(1,:)));
len=size(traces,1);

on_durations=[];
off_durations=[];
for i=1:n
    [temp_on,temp_off]=durations(traces(:,i)',len);
    on_durations=[on_durations,temp_on];
    off_durations=[off_durations,temp_off];
end

disp(n);

cd(path);
fn='on_durations.txt';
save(fn,'on_durations','-ascii');

fn='off_durations.txt';
save(fn,'off_durations','-ascii');

fn='intensity_histogram.txt';
save(fn,'countsT','-ascii');

close('all');
end







function [temp_on,temp_off]=durations(tr,len)

d_tr=conv(tr,[1,-1],'same');
d_tr(len)=0;
index1=find(d_tr==1);
index2=find(d_tr==-1);
n1=size(index2,2);
temp_on=index2-index1(1:n1);
index2=[100,index2];
n2=size(index1,2);
temp_off=index1-index2(1:n2);

end