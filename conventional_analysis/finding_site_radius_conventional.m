%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

r is the minimum separation between peaks.

Check and adjust parameters that are marked with "frank".
%}

function [good,no_good]=finding_site_radius_conventional(m,n1,n2,n3,n4,threshold,r)

good=[];

index=(m(:,1)<20.0*n1)+(m(:,1)>20.0*n2)+(m(:,2)<20.0*n3)+(m(:,2)>20.0*n4);
index=logical(index);
m(index,:)=[];

frame=unique(m(:,5));
index=m(:,5)==frame(1);
no_good=sum(index);
good=[good;m(index,1:2)];

for i = 2:length(frame)
    index=m(:,5)==frame(i);
    temp_no_good=sum(index);
    temp_good=m(index,1:2);
    for j = 1:temp_no_good
        temp=good;
        temp(:,1)=temp(:,1)-temp_good(j,1);
        temp(:,2)=temp(:,2)-temp_good(j,2);
        dist=sqrt(temp(:,1).^2+temp(:,2).^2);
        if min(dist) > 20.0*r
            no_good=no_good+1;
            good=[good;temp_good(j,:)];
        end
    end
end

for i = 1:no_good
    dist = sqrt((m(:,1)-good(i,1)).^2+(m(:,2)-good(i,2)).^2);
    temp=(dist <= 20.0*r);
    temp_x = m(temp,1);
    temp_y = m(temp,2);
    good(i,1) = mean(temp_x);
    good(i,2) = mean(temp_y);
end

index=[];
for i = 1:no_good
    dist = sqrt((m(:,1)-good(i,1)).^2+(m(:,2)-good(i,2)).^2);
    temp=(dist <= 20.0*r);
    if sum(temp) <= threshold/20.0
        index = [index;i];
    end
end
no_good=no_good-length(index);
good(index,:)=[];

good = ceil(good/20.0);

return;
end