function [id_tr,states]=analyze_localization2vbFRET(c_intensity,center_x,center_y)

x_pos=floor(center_x/20.0)+1;
y_pos=floor(center_y/20.0)+1;

FRET = cell(1,1);
FRET{1,1} = (c_intensity-min(c_intensity))/(max(c_intensity)-min(c_intensity));

[id_tr,states,K]=vbFRET_no_gui_ha(FRET,x_pos,y_pos);
id_tr=id_tr{1,K};
states=states{1,K};

end