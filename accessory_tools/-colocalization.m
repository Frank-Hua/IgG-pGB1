%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Two map files (filmx_localizations.txt) are used to calculate the colocalization percentage.

Check and adjust parameters that are marked with "frank".
%}

function colocalization

%initial input
%Get the directory and names of maps to be compared.
path=command_input('input directory','C:\\Users\\frank\\Documents\\MATLAB','s');
cd(path);
fname = command_input('index # of the reference filename ','3','s');
fname2 = command_input('index # of the data filename ','4','s');

m1 = dlmread(['film' fname '_localizations.txt']);
m2 = dlmread(['film' fname2 '_localizations.txt']);

%m2(:,1) = m2(:,1)-167; %frank
%m2(:,2) = m2(:,2)+13; %frank
%Swap columns if the movie was flipped.
temp2 = m2(:,1); %frank
m2(:,1) = m2(:,2);
m2(:,2) = temp2; %frank

%Convert unit into nm.
m1 = double(m1*20.0);
m2 = double(m2*20.0);
count0=size(m2,1);

[dist,dist_x,dist_y,avg_dist_x,avg_dist_y,count,threshold]=colocalization_core(m1,m2);
choose = command_input('perform correction?(y/n) ','y','s');

while choose=='y'
    m2(:,1) = m2(:,1)+single(avg_dist_x);
    m2(:,2) = m2(:,2)+single(avg_dist_y);
    [dist,dist_x,dist_y,avg_dist_x,avg_dist_y,count,threshold]=colocalization_core(m1,m2);
    choose = command_input('perform correction?(y/n) ','y','s');
end

%figure starts%
nbins = floor(max(dist)/20)+1;
nbins_x = floor((max(dist_x)-min(dist_x))/20)+1;
nbins_y = floor((max(dist_y)-min(dist_y))/20)+1;

figure;
input('enter-to plot ','s');
zoom on;

subplot(1,12,[1 4]);
plot(dist_x*125/180,dist_y*125/180,'x');
axis equal;

subplot(3,12,[6 10]);
hist(dist*125/180,nbins);
title('distance histogram','FontSize',9);

subplot(3,12,[18 22]);
hist(dist_x*125/180,nbins_x);
title('x distance histogram','FontSize',9);

subplot(3,12,[30 34]);
hist(dist_y*125/180,nbins_y);
title('y distance histogram','FontSize',9);

subplot(2,12,[11 12]);
text(1,6,'text window:','FontSize',15);
text(1,4,['threshold ' num2str(threshold*125/180) ' nm'],'FontSize',11);
text(1,3,['total spots ' num2str(count0)],'FontSize',11);
text(1,2,['colocalized spots ' num2str(count)],'FontSize',11);
text(1,1,'percentage of colocalization-','FontSize',11);
text(1,0,[num2str(count/size(m2,1)*100) ' %'],'FontSize',11);
temp=axis;
temp(1) = 0;
temp(2) = 8;
temp(3) = 0;
temp(4) = 7;
axis(temp);
axis off;
%figure ends%

fn=['colocalization_distance.txt'];
save(fn,'dist','-ascii');

fn=['colocalization_distance_x.txt'];
save(fn,'dist_x','-ascii');

fn=['colocalization_distance_y.txt'];
save(fn,'dist_y','-ascii');

end







%{
This function simply "tries" to calculate the offset between two sets of
    coordinates.
%}
function [dist,dist_x,dist_y,avg_dist_x,avg_dist_y,count,threshold]=colocalization_core(m1,m2)

dist = zeros(size(m2,1),1);
dist_x = zeros(size(m2,1),1);
dist_y = zeros(size(m2,1),1);
%dist_x2 = zeros(size(m2,1),1,'int32');
%dist_y2 = zeros(size(m2,1),1,'int32');
mol = zeros(size(m2,1),1,'uint16');

count = 0;
threshold = 63;

for i=1:size(m2,1)
    dist(i) = sqrt((m1(1,1)-m2(i,1))^2+(m1(1,2)-m2(i,2))^2);
    mol(i) = 1;

    for j=1:size(m1,1)
        dist2 = sqrt((m1(j,1)-m2(i,1))^2+(m1(j,2)-m2(i,2))^2);
        if dist2 < dist(i)
          dist(i) = dist2;
          mol(i) = j;
        end
    end

    if dist(i) < threshold
        count = count+1;
    end

    dist_x(i) = m1(mol(i),1)-m2(i,1);
    dist_y(i) = m1(mol(i),2)-m2(i,2);
    %dist_x2(i) = int32(round(dist_x(i)/20)*20);
    %dist_y2(i) = int32(round(dist_y(i)/20)*20);
end

nbins_x = floor((max(dist_x)-min(dist_x))/20)+1;
nbins_y = floor((max(dist_y)-min(dist_y))/20)+1;
[counts,centers] = hist(dist_x,nbins_x);
fx = fit(centers.',counts.','gauss1');
[counts,centers] = hist(dist_y,nbins_y);
fy = fit(centers.',counts.','gauss1');
avg_dist_x = fx.b1;
avg_dist_y = fy.b1;
%avg_dist_x = mode(dist_x2);
%avg_dist_y = mode(dist_y2);

disp(count);
disp(avg_dist_x);
disp(avg_dist_y);

return;
end
            
