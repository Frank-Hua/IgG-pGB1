%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Adjustable parameters are marked with "frank".
%}

function s_avg_dist=drift_correction(m)

A=dir;
[nf,~]=size(A);

%{
Go through the folder to see if "_pydrift.txt" or "_drift file.txt" exists. "_pydrift.txt" is generated by DAO-STORM. It has four columns and it is 
    in the unit of conventional pixel.
%}
g_drift = 1;
for i=1:nf
    if A(i).isdir == 0
        s=A(i).name;
        if contains(s, '_pydrift.txt')
            g_drift = 0;
            break;
        elseif contains(s, '_drift.txt')
            g_drift = -1;
            break;
        end
    end
end

%Generate the "_drift.txt" file or read existing "_*drift.txt" files.
if g_drift == 1
    a=m(:,5);
    mini = min(a);
    maxi = max(a);
    if mini ~= 0
        warning('marker does not start with frame 0');
    end
    
    temp_avg_dist=zeros(1,2);
    frame=mini;
    for i=mini+1:maxi
        index=find(a==i);
        if ~isempty(index)
            frame=[frame;i];
            dist1=m(index,1)-m(1,1);
            dist2=m(index,2)-m(1,2);
            dist=sqrt(dist1.^2+dist2.^2);
            [~,p]=min(dist);
            temp_avg_dist(end+1,1)=dist1(p);
            temp_avg_dist(end,2)=dist2(p);
        end
    end
    avg_dist(:,1)=pchip(frame,temp_avg_dist(:,1),[mini:maxi]');
    avg_dist(:,2)=pchip(frame,temp_avg_dist(:,2),[mini:maxi]');
elseif g_drift == 0
    temp_avg_dist=dlmread(s);
    avg_dist(:,1)=180.0*temp_avg_dist(:,2);
    avg_dist(:,2)=180.0*temp_avg_dist(:,3);

elseif g_drift == -1
    temp_avg_dist=dlmread(s);
    avg_dist(:,1)=temp_avg_dist(:,1);
    avg_dist(:,2)=temp_avg_dist(:,2);
end

hdl = figure;
drift_plot(hdl,avg_dist,'b','b');
s_avg_dist(:,1)=smooth(avg_dist(:,1),21,'sgolay');
s_avg_dist(:,2)=smooth(avg_dist(:,2),21,'sgolay');
drift_plot(hdl,s_avg_dist,'k','k');
input('enter-to continue ','s');
close(hdl);

return;
end







function drift_plot(hdl,avg_dist,c1,c2)

figure(hdl);
subplot(2,1,1);
plot(avg_dist(:,1)/180.0,c1);
title ('x drift (conventional pixels)');
hold on;
zoom on;
subplot(2,1,2);
plot(avg_dist(:,2)/180.0,c2);
title ('y drift (conventional pixels)');
hold on;
zoom on;

return;
end