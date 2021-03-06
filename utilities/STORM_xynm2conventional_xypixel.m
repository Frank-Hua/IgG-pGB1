%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Adjustable parameters are marked with "frank".
%}

function [fx_pos,fy_pos]=STORM_xynm2conventional_xypixel(center_x,center_y,situ)

switch situ
    case 1
        %Use these if conventional image is 256*256 pixel^2 and xy flipped.
        fx_pos = center_y/180.0;
        fy_pos = center_x/180.0-1.0;
    case 2
        %Use these if conventional image is 512*512 pixel^2.
        fx_pos = center_x/180.0+127.0;
        fy_pos = center_y/180.0+128.0;
    case 3
        %Use these if conventional image is 512*512 pixel^2 and xy flipped.
        fx_pos = center_y/180.0+128.0;
        fy_pos = center_x/180.0+127.0;
    otherwise
        %Use these if conventional image is 256*256 pixel^2.
        fx_pos = center_x/180.0-1.0;
        fy_pos = center_y/180.0;
end