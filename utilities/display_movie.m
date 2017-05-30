%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function I=display_movie(frame,tr2,len,center_x,center_y,s_avg_dist,situ)

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

img_size = 17; %frank
s_avg_dist = flip_drift_correction(s_avg_dist,mod(situ,2));
s_avg_dist = round(s_avg_dist/180.0);

%{
center_xy is a floating number recording xy coordinate in the unit of nm.
fxy_pos is a floating number recording xy coordinate in the unit of
    diffration-limited pixel number.
xy_pos is an integer number recording xy coordinate in the unit of
    diffration-limited pixel number.
%}
[fx_pos,fy_pos] = STORM_xynm2conventional_xypixel(center_x,center_y,situ);
x_pos = floor(fx_pos)+1;
y_pos = floor(fy_pos)+1;

index=find(tr2);
count=length(index);
I = zeros(img_size,img_size,count,'uint16');

for n=1:count
    t=index(n);
    %avg_dist(t,i) is the drift correction at frame t, with i=1 for x and i=2 for y.
    I(:,:,n) = frame(x_pos+s_avg_dist(t,1)-(img_size-1)/2:x_pos+s_avg_dist(t,1)+(img_size-1)/2,y_pos+s_avg_dist(t,2)-(img_size-1)/2:y_pos+s_avg_dist(t,2)+(img_size-1)/2,t);
end

I = mat2gray(I);

return;
end