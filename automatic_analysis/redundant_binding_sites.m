%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

This is for skipping binding sites that have been analyzed.

Check and adjust parameters that are marked with "frank".
%}

function [good,no_good] = redundant_binding_sites(path,good,no_good)

cd([path '\molecules']);
A=dir;
[nf,~]=size(A);

for i=1:nf
    s1=A(i).name;
    if A(i).isdir == 0 && strcmp(s1(end-2:end), 'txt')
        %Each file is named as something_x-y-#.txt
        p1 = strfind(s1,'_');
        p1 = [p1 strfind(s1,'-')];
        if ~isempty(p1)
            x_pos1=str2double(s1(p1(1)+1:p1(2)-1));
            y_pos1=str2double(s1(p1(2)+1:p1(3)-1));

            dist = sqrt((good(:,1)-x_pos1).^2+(good(:,2)-y_pos1).^2);
            [dist0,index] = min(dist);

            if dist0 < 3
                good(index,:)=[];
                no_good = no_good-1;
            end
        end
    end
end

cd(path);

return;
end