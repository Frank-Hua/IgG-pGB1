%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

total records background intensity, KM, kon and koff of each binding site.
intensity2 records intensities when protein G is bound to a binding site.
tr records the digital trace of each binding site.
value determines if the analyzed values above get accepted.

For some reason, the pixel size for diffration-limited image is set as 180.0 nm. And the "pixel size" for STORM image is 20.0 nm.

This version is for fast analysis by skipping manual binding site picking.

This version is for conventional analysis.

Adjustable parameters are marked with "frank".
%}

function [total,intensity2,value]=analyze_site_conventional(n,m,center_x,center_y,frame,len,s_avg_dist,baseline,hdl,hdl2,situ)

total=zeros(4,1);
intensity2=0;
value=0;

%{
This is to find localization events within a certain radius of center_xy.

    flag stores the number of found localizations.
    local_xy stores the xy coordinate of found localizations.
%}
dist = sqrt((m(:,1)-center_x).^2+(m(:,2)-center_y).^2);
index=(dist <= 450);
flag=sum(index);
flag_num = num2str(flag);
local_x = m(index,1);
local_y = m(index,2);

answer = '';
%280.0 nm is the radius for defining one binding site.
rad = 280.0; %frank
rad_num = num2str(125.0/180.0*rad); %frank; correction for pixel size

%Frame number is used as timeunit.
timeunit = 1; %frank
time = (1:len)*timeunit;
clicktime2 = [];
clicktime3 = [];
thres = 10000;

while ~strcmp(answer,'done')
    %smm_intensity function analyzes the intensity of a binding site.
    intensity=smm_intensity(frame,len,center_x,center_y,s_avg_dist,situ);
    nbins = floor((max(intensity))/50);
    [counts,centers] = hist(intensity,nbins);
    bsl0 = mean(intensity(1:20));
    try
        f = fit(centers.',counts.','gauss1');
        bsl  = f.b1;
    catch
        bsl = bsl0;
    end
    
    %Correct intensity baseline.
    baseline(:,1)=baseline(:,1)/mean(baseline(:,1))*bsl;
    c_intensity = intensity./baseline(:,1)*bsl;
    
    %This is to find frames with intensity above threshold.
    f_intensity=ordfilt2(c_intensity,3,[1;1;1]);

    tr=f_intensity>thres;
    tr=tr*1.0;
    tr=imerode(tr,strel('line',2,90));
    
    %{
    Fill up 2-frame gaps in tr that are caused by incomplete STROM localization detection.
    %}
    tr=imdilate(tr,strel('line',2,90));
    tr=imerode(tr,strel('line',2,90));
    
    %{
    clicktime2 and clicktime3 are used for correcting tr. And clicktime2 is always overwritten by clicktime3.
    %}
    if length(clicktime2)>=2
        for t=1:2:(length(clicktime2)-1)
            tr(floor(clicktime2(t)/timeunit):floor(clicktime2(t+1)/timeunit)) = 1;
        end
    end
    if length(clicktime3)>=2
        for t=1:2:(length(clicktime3)-1)
            tr(floor(clicktime3(t)/timeunit):floor(clicktime3(t+1)/timeunit)) = 0;
        end
    end
    
    %inCircle_num stores the number of localizations within a binding site.
    aboveThres_num = num2str(nnz(tr));
    
    %Generate a circle around a binding site.
    circle = zeros(180,2);
    circle(:,1) = rad*cos(2*(1:180)/180*pi)+center_x;
    circle(:,2) = rad*sin(2*(1:180)/180*pi)+center_y;
    
    %Just plot segments of tr when protein G was bound.
    t_tr = imdilate(tr,strel('line',50,90));
    index = find(t_tr);
    
    %Generate indicator lines based on index.
    if ~isempty(index)
        [line_x,line_y]=line_generator(index);
    else
        line_x=[0;0];
        line_y=[0;1];
    end
    
    %figure starts%
    figure(hdl);
    subplot(3,10,[1 6],'replace');
    plot(time(1:floor(len/2)),c_intensity(1:floor(len/2)),'b');
    temp=plot_formatter('time trace 1st half (intensity+presence)',1,1,1,1,'off','off',0,max(c_intensity)*1.1);
    plot(time(1:floor(len/2)),0.25*temp(4)*tr(1:floor(len/2)),'r');

    subplot(3,10,[11 16],'replace');
    plot(time(floor(len/2)+1:len),c_intensity(floor(len/2)+1:len),'b');
    temp=plot_formatter('time trace 2nd half (intensity+presence)',1,1,1,1,'off','off',0,max(c_intensity)*1.1);
    plot(time(floor(len/2)+1:len),0.25*temp(4)*tr(floor(len/2)+1:len),'r');

    subplot(3,10,[21 26],'replace');
    plot(1:length(index),c_intensity(index),'b');
    temp=plot_formatter('time trace (protein G-bound portion) (intensity+presence)',1,1,1,1,'off','off',0,max(c_intensity)*1.1);
    plot(1:length(index),0.25*temp(4)*tr(index),'r');
    line(line_x,max(c_intensity)*1.1*line_y,'LineStyle','--');
    
    %Baseline is corrected, so (0,bsl) is plotted.
    subplot(3,10,17,'replace');
    plot(0,bsl,'bo');
    plot_formatter('',0,0,0,0,'off','off',0,max(c_intensity)*1.1);
    
    subplot(2,10,[8 10],'replace');
    text(1,5,['text window of binding site ' num2str(n) ':'],'FontSize',15);
    text(1,4,['total number of localizations: ' flag_num],'FontSize',11);
    text(1,3,['the radius is ' rad_num ' nm'],'FontSize',11);
    text(1,2,['number of detected frames:' aboveThres_num],'FontSize',11);
    plot_formatter('',0,0,0,0,0,10,0,6);

    subplot(2,10,[18 20],'replace');
    plot(local_x,local_y,'+k',circle(:,1),circle(:,2));
    axis equal;
    axis square;
    temp=axis;
    temp(1) = center_x-500; temp(2) = center_x+500; temp(3) = center_y-500; temp(4) = center_y+500;
    axis(temp);
    %figure ends%
    
    X = [];
    Y = [];
    display_menu2;
    answer = input('option: ','s');
    
    %Option 'c' is to adjust the center of a binding site.
    if answer=='c'
        [X,Y] = ginput(1);
        center_x = X;
        center_y = Y;
    end
    
    %Option 'd' is to calculate the distance between any two points.
    if answer=='d'
        [X,Y] = ginput(2);
        dist2 = sqrt((X(1)-X(2))^2+(Y(1)-Y(2))^2);
        dist2_num = num2str(125.0/180.0*dist2); %correction for pixel size
        
        %figure starts%
        subplot(2,10,[8 10]);
        text(1,1,['the distance is ' dist2_num ' nm'],'FontSize',11);
        %figure ends%
        
        input('enter-to continue ','s');
    end
    
    %Option 't' is to change intensity threshold.
    if answer=='t'
        [~,thres] = ginput(1);
    end
    
    %Option 'p' is to display the raw intensity and baseline.
    if answer=='p'
        %figure starts%
        figure(hdl);
        subplot(3,10,[1 6]);
        plot(time(1:floor(len/2)),intensity(1:floor(len/2)),'k');
        plot_formatter('time trace 1st half (intensity+presence)',1,1,1,0,'off','off',0,max(c_intensity)*1.1);
        plot(time(1:floor(len/2)),baseline(1:floor(len/2),1),'g');

        subplot(3,10,[11 16]);
        plot(time(floor(len/2)+1:len),intensity(floor(len/2)+1:len),'k');
        plot_formatter('time trace 2nd half (intensity+presence)',1,1,1,0,'off','off',0,max(c_intensity)*1.1);
        plot(time(floor(len/2)+1:len),baseline(floor(len/2)+1:len,1),'g');

        subplot(3,10,17);
        plot(min(counts)+20,bsl,'k+',counts,centers,'bo');
        plot_formatter('',0,0,0,0,'off','off',0,max(c_intensity)*1.1);
        %figure ends%

        input('enter-to continue ','s');
    end
    
    %Option 'f' is to modify tr.
    if answer=='f'
        disp('left click to add frames');
        disp('right click to delete frames');
        
        %figure starts%
        clf(hdl2,'reset')
        figure(hdl2);
        plot(1:length(index),c_intensity(index),'b');
        temp=plot_formatter('time trace (protein G-bound portion) (intensity+presence)',1,1,1,1,'off','off',0,max(c_intensity)*1.1);
        plot(1:length(index),0.25*temp(4)*tr(index),'r');
        line(line_x,max(c_intensity)*1.1*line_y,'LineStyle','--');
        %figure ends%
        
        [clicktime2,clicktime3]=tr_modifier(index,clicktime2,clicktime3);
    end
    
    %Option 'q' is to analyze the current binding site and exit.
    if answer=='q'
        intensity2=nonzeros(c_intensity.*tr);

        d_tr=conv(tr,[1,-1],'same');
        d_tr(len)=0;
        index=find(d_tr==1);
        total(1)=bsl-bsl0;
        total(2)=index(end)/nnz(tr); %frank KM = total(2)*conc
        total(3)=size(index,1)/index(end); %frank kon=total(3)/(exposure_time*conc)
        total(4)=total(2)*total(3); %frank koff=total(4)/exposure_time
        
        value=-1;
        answer='done';
    end
    
    %Option 'g' is to navigate to a different binding site by its number.
    if answer=='g'
        value=input('which binding site: ');
        answer='done';
    end
    
    %Option 'exit' is to terminate analysis.
    if strcmp(answer,'exit')
        value=-2;
        answer='done';
    end
    
end

return;
end







function display_menu2

disp('======================================================================');
disp('& analysis menu options &');
disp('center-(c), distance-(d), threshold-(t)');
disp('check the raw intensity and baseline-(p)');
disp('modifiy the time trace (presence)-(f), analyze-(q)');
disp('navigate to a different binding site-(g)');
disp('stop analyzing the current binding site-(done)');
disp('stop analyzing-(exit)');

return;
end