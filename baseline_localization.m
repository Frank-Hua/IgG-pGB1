%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

baseline records intensity baseline and it is generated from a blank area
    where hardly any protein G bound.

For some reason, the pixel size for diffration-limited image is set as 180.0 nm.
The "pixel size" for STORM image is 20.0 nm.

Check and adjust parameters that are marked with "frank".
%}

function baseline = baseline_localization(n,m,center_x,center_y,frame,len,s_avg_dist,hdl,situ)

baseline = [];
%{
This is to find localization events within a certain radius of center_xy.

    flag stores the number of found localizations.
    local_xy stores the xy coordinate of found localizations.
    frame_num stores the frame number of a frame in which a localziation
        is detected.
%}
dist = sqrt((m(:,1)-center_x).^2+(m(:,2)-center_y).^2);
index=(dist <= 450);
flag=sum(index);
flag_num = num2str(flag);
local_x = m(index,1);
local_y = m(index,2);
frame_num = m(index,5)+1;

answer = '';
%70.0 nm is the radius for defining one binding site.
rad = 70.0; %frank
rad_num = num2str(125.0/180.0*rad); %correction for pixel size

%Frame number is used as timeunit.
timeunit = 1; %frank
time = (1:len)*timeunit;

while ~strcmp(answer,'done')
    %smm_intensity function analyzes the intensity of a binding site.
    intensity=smm_intensity(frame,len,center_x,center_y,s_avg_dist,situ);
    nbins = floor((max(intensity))/50);
    [counts,centers] = hist(intensity,nbins);
    f = fit(centers.',counts.','gauss1');
    bsl0 = mean(intensity(1:20));
    bsl  = f.b1;
    
    %This is to find localization events from one binding site and generate tr.
    tr = zeros(len,1);
    dist = sqrt((local_x-center_x).^2+(local_y-center_y).^2);
    index=(dist <= rad);
    tr(frame_num(index))=1;
    
    %{
    Fill up 2-frame gaps in tr that are caused by incomplete STROM
        localization detection.
    %}
    tr=imdilate(tr,strel('line',2,90));
    tr=imerode(tr,strel('line',2,90));
    
    %inCircle_num stores the number of localizations within a binding site.
    inCircle_num = num2str(nnz(tr));
    
    %Correct intensity baseline.
    [c_intensity,t_baseline] = baseline_correction(intensity,tr,bsl0,bsl);
    
    %Generate a circle around a binding site.
    circle = zeros(180,2);
    circle(:,1) = rad*cos(2*(1:180)/180*pi)+center_x;
    circle(:,2) = rad*sin(2*(1:180)/180*pi)+center_y;
    
    %figure starts%
    figure(hdl);
    
    subplot(2,10,[1 6],'replace');
    plot(time(1:floor(len/2)),c_intensity(1:floor(len/2)),'b');
    temp=plot_formatter('time trace 1st half (intensity+presence)',1,1,1,1,'off','off',0,max(c_intensity)*1.1);
    plot(time(1:floor(len/2)),0.25*temp(4)*tr(1:floor(len/2)),'r');
    plot(time(1:floor(len/2)),intensity(1:floor(len/2)),'k');
    plot(time(1:floor(len/2)),t_baseline(1:floor(len/2)),'g');
    
    subplot(2,10,[11 16],'replace');
    plot(time(floor(len/2)+1:len),c_intensity(floor(len/2)+1:len),'b');
    temp=plot_formatter('time trace 2nd half (intensity+presence)',1,1,1,1,'off','off',0,max(c_intensity)*1.1);
    plot(time(floor(len/2)+1:len),0.25*temp(4)*tr(floor(len/2)+1:len),'r');
    plot(time(floor(len/2)+1:len),intensity(floor(len/2)+1:len),'k');
    plot(time(floor(len/2)+1:len),t_baseline(floor(len/2)+1:len),'g');
        
    subplot(2,10,17,'replace');
    plot(min(counts)+20,bsl,'k+',counts,centers,'bo');
    plot_formatter('',0,0,0,0,'off','off',0,max(c_intensity)*1.1);
    
    subplot(2,10,[7 10],'replace');
    text(1,5,['text window of localization ' num2str(n) ':'],'FontSize',15);
    text(1,4,['total number of localizations: ' flag_num],'FontSize',11);
    text(1,3,['the radius is ' rad_num ' nm'],'FontSize',11);
    text(1,2,['number of localizations in the circle:' inCircle_num],'FontSize',11);
    plot_formatter('',0,0,0,0,0,10,0,6);

    subplot(2,10,[18 20],'replace');
    plot(local_x,local_y,'+k',circle(:,1),circle(:,2));
    axis equal;
    axis square;
    temp=axis;
    temp(1) = center_x-500; temp(2) = center_x+500; temp(3) = center_y-500; temp(4) = center_y+500;
    axis(temp);
    %figure ends%
    
    X = zeros;
    Y = zeros;
    display_menu2;    
    answer = input('option: ','s');
    
    %Option 'c' is to adjust the center of a binding site.
    if answer=='c'
        [X,Y] = ginput(1);
        center_x = X;
        center_y = Y;
    end
    
    %Option 'r' is to adjust the radius that defines a binding site.
    if answer=='r'
        [X,Y] = ginput(1);
        rad = sqrt((X-center_x)^2+(Y-center_y)^2);
        rad_num = num2str(125.0/180.0*rad); %correction for pixel size
    end
    
    %Option 'd' is to calculate the distance between any two points.
    if answer=='d'
        [X,Y] = ginput(2);
        dist2 = sqrt((X(1)-X(2))^2+(Y(1)-Y(2))^2);
        dist2_num = num2str(125.0/180.0*dist2); %correction for pixel size
        
        %figure starts%
        subplot(2,10,[7 10]);
        text(1,1,['the distance is ' dist2_num ' nm'],'FontSize',11);
        %figure ends%
        
        input('enter-to continue ','s');
    end
    
    %Option 'b' is to confirm the current baseline and exit.
    if answer=='v'
        baseline=t_baseline;
        answer='done';
    end
end

return;
end







function display_menu2

disp('======================================================================');
disp('& baseline correction menu options &');
disp('center-(c), radius-(r), distance-(d), confirm baseline-(v)');
disp('stop generating baseline from the current localization-(done)');

return;
end