%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

This is for skipping redundant binding sites.

Check and adjust parameters that are marked with "frank".
%}

function [good,no_good] = redundant_binding_site(path,good,no_good)

cd([path '\molecules']);
A=dir;
[nf,dum]=size(A);

for i=1:nf
    s1=A(i).name;
    if A(i).isdir == 0 && strcmp(s1(end-2:end), 'txt')
        p1 = strfind(s1,'_');
        p1 = [p1 strfind(s1,'-')];
        x_pos1=str2double(s1(p1(1)+1:p1(2)-1));
        y_pos1=str2double(s1(p1(2)+1:p1(3)-1));
        dist0=10000;
        for j=1:no_good
            x_pos2=good(j,1);
            y_pos2=good(j,2);
            dist = sqrt((x_pos1-x_pos2)^2+(y_pos1-y_pos2)^2);
            if dist < dist0
                dist0 = dist;
                index = j;
            end
        end
        if dist0 < 3
            good(index,:)=[];
            no_good = no_good-1;
        end
    end
end

cd(path);

return;
end