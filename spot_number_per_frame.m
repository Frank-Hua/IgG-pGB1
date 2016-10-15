%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function spot_number_per_frame(fname,m)

mini=min(m(:,5));
maxi=max(m(:,5));
count=zeros((maxi-mini+1),1);
delete(['film' fname '_snpf.txt']);
fid = fopen( ['film' fname '_snpf.txt'], 'a' );
j=1;
for i=mini:maxi
    while m(j,5) == i
        count(i-mini+1)=count(i-mini+1)+1;
        j=j+1;
        if j > size(m,1)
            break
        end
    end
    fprintf(fid,'%f\t%f\n',i,count(i-mini+1));
end
fclose(fid);

hdl = figure;
snpf_plot(hdl,count);
input('enter-to continue ','s');
close(hdl);

return;
end







function snpf_plot(hdl,count)

figure(hdl);
plot(count);
title ('spot number per frame');
hold on;
zoom on;
count2=smooth(count,51,'sgolay');
plot(count2,'k');

return;
end