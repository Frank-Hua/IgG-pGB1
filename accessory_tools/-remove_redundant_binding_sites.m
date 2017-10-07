%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

This is for scanning analyzed binding sites and getting rid of ones that 
    appeared more than once.

Check and adjust parameters that are marked with "frank".
%}

function remove_redundant_binding_sites

%initial input
%Get the directory that needs to be scanned.
path=command_input('input directory','C:\\Users\\frank\\Documents\\MATLAB','s');

for t=1:3
    %{
    Under each directory, there are three folders to be scanned. We will
        go through each of them.
    %}
    switch t
        case 1
            cd([path '\molecules']);
        case 2
            cd([path '\intensities']);
        case 3
            cd([path '\traces']);
    end

    A=dir;
    [nf,dum]=size(A);
    
    %{
    Each possible pair of two files are compared to determine if redundancy 
        exists.
    %}
    for i=1:nf-1
        for j=i+1:nf
            s1=A(i).name;
            s2=A(j).name;
            if A(i).isdir == 0 && A(j).isdir == 0 && strcmp(s1(end-2:end), 'txt') && strcmp(s2(end-2:end), 'txt')
                %Each file is named as something_x-y-#.txt
                p1 = strfind(s1,'_');
                p1 = [p1 strfind(s1,'-')];
                p2 = strfind(s2,'_');
                p2 = [p2 strfind(s2,'-')];
                x_pos1=str2double(s1(p1(1)+1:p1(2)-1));
                y_pos1=str2double(s1(p1(2)+1:p1(3)-1));
                x_pos2=str2double(s2(p2(1)+1:p2(2)-1));
                y_pos2=str2double(s2(p2(2)+1:p2(3)-1));
                dist = sqrt((x_pos1-x_pos2)^2+(y_pos1-y_pos2)^2);
                %Compare the date created and delete the older one.
                if dist < 3
                    if A(i).datenum <= A(j).datenum
                        disp(A(i).name);
                        delete(A(i).name);
                    else
                        disp(A(j).name);
                        delete(A(j).name);
                    end
                end
            end
        end
    end

    cd(path);
end
