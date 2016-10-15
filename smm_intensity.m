%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function intensity=smm_intensity(frame,len,center_x,center_y,s_avg_dist,situ)

%{
switch situ
	case 0
        %Use these if conventional image is 256*256 pixel^2.
    case 1
        %Use these if conventional image is 256*256 pixel^2 and xy flipped.
    case 2
        %Use these if conventional image is 512*512 pixel^2.
    case 3
        %Use these if conventional image is 512*512 pixel^2 and xy flipped.
%}

intensity = zeros(len,1);
s_avg_dist=flip_drift_correction(s_avg_dist,mod(situ,2));
s_avg_dist = round(s_avg_dist/180.0);

%{
center_xy is a floating number recording xy coordinate in the unit of nm.
fxy_pos is a floating number recording xy coordinate in the unit of
    diffration-limited pixel number.
xy_pos is an integer number recording xy coordinate in the unit of
    diffration-limited pixel number.
%}
[fx_pos,fy_pos]=STORM_xynm2conventional_xypixel(center_x,center_y,situ);
x_pos = floor(fx_pos)+1;
y_pos = floor(fy_pos)+1;

g_peaks = zeros(3,3,7,7);
for k = 1:3
    for l = 1:3
        offx = 0.5*(k-2);
        offy = 0.5*(l-2);
        for i = 1:7
            for j = 1:7
                dist = 0.4*((i-4-offx)^2+(j-4-offy)^2); %frank
                g_peaks(k,l,i,j) = exp(-dist);
            end
        end
    end
end
k = round(2*(fx_pos-x_pos+0.5))+2;
l = round(2*(fy_pos-y_pos+0.5))+2;
%dum is used for gaussian weighting.
dum = reshape(g_peaks(k,l,:,:),7,7);

ave_arr=zeros(7,7,'uint16');

for t=1:10
    ave_arr = ave_arr + frame(x_pos-3:x_pos+3,y_pos-3:y_pos+3,t);
end
ave_arr = ave_arr/10;
background = min(min(ave_arr)); %frank

for t=1:len
    %dum2 = sum(frame(x_pos+s_avg_dist(t,2)-3:x_pos+s_avg_dist(t,2)+3,y_pos+s_avg_dist(t,1)-3:y_pos+s_avg_dist(t,1)+3,t));
    %intensity(t) = sum(dum2)-49*background;
    
    %avg_dist(t,i) is the drift correction at time t, with i=1 for x and i=2 for y.
    dum1 = double(frame(x_pos+s_avg_dist(t,1)-3:x_pos+s_avg_dist(t,1)+3,y_pos+s_avg_dist(t,2)-3:y_pos+s_avg_dist(t,2)+3,t)-background);
    dum2 = sum(dum.*dum1);
    intensity(t) = sum(dum2);
end

return;
end