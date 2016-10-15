%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function drift_correction_control(fname,m,m0)

s_avg_dist=drift_correction(m0);
delete(['film' fname '_drift.txt']);
fid = fopen(['film' fname '_drift.txt'], 'a' );
for j=1:size(s_avg_dist,1)
    fprintf(fid,'%f\t%f\n',s_avg_dist(j,:));
end
fclose(fid);

m=drift_apply(m,s_avg_dist);
delete(['film' fname '_drift_corrected.txt']);
fid = fopen(['film' fname '_drift_corrected.txt'],'a');
for j=1:size(m,1)
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\n',m(j,:));
end
fclose(fid);

count_size = 2295; %frank
count=histogram_2D(m,count_size);
delete(['film' fname '_drift_corrected_histogram.tif']);
imwrite(uint16(count),['film' fname '_drift_corrected_histogram.tif']);

return;
end







function mo=drift_apply(mi,s_avg_dist)

%mini = min(mi(:,5));
j=1;
while j <= size(mi,1)
    %mi(j,1) = mi(j,1)-s_avg_dist(mi(j,5)-mini+1,1);
    %mi(j,2) = mi(j,2)-s_avg_dist(mi(j,5)-mini+1,2);
    mi(j,1) = mi(j,1)-s_avg_dist(mi(j,5)+1,1);
    mi(j,2) = mi(j,2)-s_avg_dist(mi(j,5)+1,2);
    j = j+1;
end
mo=mi;

return;
end

function count=histogram_2D(m,count_size)

count=zeros(count_size,count_size);
for k=1:size(m,1)
    i = floor(m(k,1)/20.0)+1;
    j = floor(m(k,2)/20.0)+1;
    if i >= 1 && i <= count_size && j >= 1 && j <= count_size
        count(i,j) = count(i,j)+1;
    end
end
count = uint32(count);
count = 20*count;

return;
end