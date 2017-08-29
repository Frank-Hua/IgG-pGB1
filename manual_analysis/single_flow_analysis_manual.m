%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Output files are "molecule_xx-yy-#", "intensity_xx-yy-#" and "trace_xx-yy-#".

    "molecule_xx-yy-#" files record background intensity, KM, kon and koff.
    "intensity_xx-yy-#" files record intensities when protein G is bound.
    "trace_xx-yy-#" files record the digital trace.

For some reason, the pixel size for diffration-limited image is set as 180.0 nm. And the "pixel size" for STORM image is 20.0 nm.

Adjustable parameters are marked with "frank".
%}

function single_flow_analysis_manual

%% initial input
%command_input function doesn't do much, but is good to have.
addpath('D:\Hua\02documents\08MatLab scripts\protein G\IgG-pGB1');
addpath('D:\Hua\02documents\08MatLab scripts\protein G\IgG-pGB1\utilities');

path=command_input('input directory:','C:\\Users\\frank\\Documents\\MATLAB','s');
cd(path);
fname = command_input('input file index #:','1','s');
%Choose the size of the frame and whether to transpose the frame.
situ = input('movie option: [(0)-256, (1)-flipped 256, (2)-512, (3)-flipped 512] ');
if isempty(situ)
    situ=0;
end


%% read smm files and other image info
%.smm file (turned into frame) provides intensity information of selected binding sites.
fid = fopen(['film' fname '.smm'],'r');
if ( fid == -1 )
	display( 'file open failed' );
	return;
end
fileinfo = dir(['film' fname '.smm']);
film_x = fread( fid, 1, 'uint16' );
film_y = fread( fid, 1, 'uint16' );
bpp = fread( fid, 1, 'uint8' );
background = fread( fid, 1, 'uint32' );
data_scaler = fread( fid, 1, 'uint32' );
framecycle = fread( fid, 1, 'float32' );
headersize = 2+2+1+4+4+4;
len = uint32( ( fileinfo.bytes - headersize ) * 1.0 / bpp / film_x / film_y );

frame = zeros(film_x,film_y,len,'uint16');
for t=1:len
    frame(:,:,t) = fread(fid,[film_x,film_y],'uint16');
end
fclose(fid);


%% perform drift correction
%{
"_drift_corrected.txt" files can be directly generated from IDL STORM codes.
    In such a case, "_drift_corrected.txt" is the same as "_raw.txt", and you can feed drift_correction_control function with a mock "_marker.txt" 
    file. And "_drift.txt" files need to be generated separately.

"_marker.txt" files must contain the first and last frames.
%}
answer = input('skip performing drift correction: [yes-(enter) or no-(no)] ','s');
if strcmp(answer,'no')
    m = dlmread(['film' fname '_raw.txt']);
    m(:,5) = uint16(m(:,5));
    m0 = dlmread(['film' fname '_marker.txt']);
    m0(:,5) = uint16(m0(:,5));
    drift_correction_control(fname,m,m0);
end

%{
"_drift_corrected_histogram.tif" is basically a STORM image of all binding sites. Its pixel size is 9 times smaller than that of a diffration-limited 
    image.
"_drift.txt" is a drift file. Drifts in x and y are recorded in the unit of nm.
%}
m2 = dlmread(['film' fname '_drift_corrected.txt']);
m2(:,5) = uint16(m2(:,5));
m3 = imread(['film' fname '_drift_corrected_histogram.tif'],'TIFF');
m4 = dlmread(['film' fname '_drift.txt']);


%% plot spot number per frame
%This can be used to calculate photobleaching rates.
answer2 = input('skip calculating spot number per frame: [yes-(enter) or no-(no)] ','s');
if strcmp(answer2,'no')
    spot_number_per_frame(fname,m2);    
end


%% start manual selection and analysis
baseline = ones(len,1);

mkdir('molecules');
mkdir('intensities');
mkdir('traces');

answer3 = command_input('specify binding site finding parameter-circle radius in sr-pixels:','5','');

answer4 = command_input('specify binding site finding parameter-minimum localization event number:','100','');
%For some reason, this value needs to be multiplied by 20.
answer4 = answer4*20.0;

n = 0;
hdl  = figure;
hdl2 = figure;
hdl3 = figure;
x = 8;
y = 8;
answer = '';

while ~strcmp(answer,'exit')
    %Display binding sites in the left window.
    %{
    The whole image is divided in to 225=15*15 windows and each window is displayed individually. Each window has a size of 153 STORM pixels or 17 
        diffraction-limited pixels.
    %}
    low_x = x*3060;
    low_y = y*3060;
    high_x = (x+1)*3060;
    high_y = (y+1)*3060;
    x_pos=floor(high_x/20.0)-76;
    y_pos=floor(high_y/20.0)-76;
    %{
    6 and 148 set the boundary for spot finding. It avoids spots close to the edges.
    %}
    [good,no_good] = finding_site_radius(m3(x_pos-76:x_pos+76,y_pos-76:y_pos+76),6,148,6,148,answer4,answer3);
    
    %figure starts%
    figure(hdl);
    subplot(1,12,[1 5]);
    hold on;
    image(m3(x_pos-76:x_pos+76,y_pos-76:y_pos+76)');
    %colorbar;
    plot(good(:,1),good(:,2),'+b');
    title(['number of binding sites: ' num2str(no_good)]);
    axis equal;
    temp=axis;
    temp(1) = -4;
    temp(2) = 158;
    temp(3) = -4;
    temp(4) = 158;
    axis(temp);
    zoom on;
    hold off;
    
    subplot(2,12,[6 7]);
    plot(x,y,'ks');
    for j=0:15
        line([j-0.5 j-0.5],[-0.5 14.5]);
        line([-0.5 14.5],[j-0.5 j-0.5]);
    end
    title(['overview: (' num2str(x) ',' num2str(y) ')']);
    temp=axis;
    temp(1) = -0.5;
    temp(2) = 14.5;
    temp(3) = -0.5;
    temp(4) = 14.5;
    axis(temp);
    axis equal;
    axis off;
    %figure ends%
    
    
    %input options
    display_menu;
    answer = input('option: ','s');
    
    %Options 'a', 'd', 's' and 'w' are to navigate display window.
    if answer == 'a'
        x = x-1;
        if x < 0
            x = 0;
        end
    end
    if answer == 'd'
        x = x+1;
        if x > 14
            x = 14;
        end
    end
    if answer == 's'
        y = y-1;
        if y < 0
            y = 0;
        end
    end
    if answer == 'w'
        y = y+1;
        if y > 14
            y = 14;
        end
    end
    
    %Option 'b' is to correct trace baseline caused by focus drift.
    if answer == 'b'
        [center_x,center_y]=ginput(1);
        %Convert center_xy into values in the unit of nm.
        center_x=low_x+20.0*center_x;
        center_y=low_y+20.0*center_y;
        
        t_baseline = baseline_site(n,m2,center_x,center_y,frame,len,m4,hdl3,situ);
        baseline = [baseline t_baseline];
        if size(baseline,2) > 1
            baseline(:,1) = mean(baseline(:,2:end),2);
        end
    end
    
    %Option 'k' is to zoom in to a binding site.
    if answer == 'k'
        [center_x,center_y]=ginput(1);
        %Convert center_xy into values in the unit of nm.
        center_x=low_x+20.0*center_x;
        center_y=low_y+20.0*center_y;
        [X,Y,~]=peaks(11);
        
        %Get xy_pos in the unit of STORM pixel number.
        x_pos=floor(center_x/20.0)+1;
        y_pos=floor(center_y/20.0)+1;
        Z=double(m3(x_pos-5:x_pos+5,y_pos-5:y_pos+5));
        
        %figure starts%
        figure(hdl);
        subplot(1,12,[8 12]);
        surf(X,Y,Z);
        rotate3d on;
        answer = input('precede to analysis: yes-(enter) or no-(no) ','s');
        %figure ends%
        
        if ~strcmp(answer,'no')
            [total,intensity2,tr,value]=analyze_site_manual(n,m2,center_x,center_y,frame,len,m4,baseline,hdl,hdl2,situ);
            if value == -1
                save_molecule(total,x_pos,y_pos,n);
                save_intensity2(intensity2,x_pos,y_pos,n);
                save_tr(tr,x_pos,y_pos,n);
            end
        end
    end
    clf(hdl,'reset');
end

close('all');
combine_molecules;
combine_intensities
combine_traces;

delete(['film' fname '_baseline.txt']);
fid = fopen(['film' fname '_baseline.txt'], 'a' );
for j=1:size(baseline(:,1),1)
    fprintf(fid,'%f\n',baseline(j,1));
end
fclose(fid);

end







function display_menu

disp('======================================================================');
disp('& main menu options &');
disp('to move the window: left-(a), right-(d), down-(s) or up-(w)');
disp('to select a site to correct intensity baseline-(b)');
disp('to select a site to analyze-(k)');
disp('to stop analyzing-(exit)');

return;
end

function save_molecule(total,x_pos,y_pos,i)

total;
cd('molecules');
fn=['molecule_' num2str(x_pos) '-' num2str(y_pos) '-' num2str(i) '.txt'];
save(fn,'total','-ascii');
cd('..');

return;
end

function save_intensity2(intensity2,x_pos,y_pos,i)

intensity2;
cd('intensities');
fn=['intensity_' num2str(x_pos) '-' num2str(y_pos) '-' num2str(i) '.txt'];
save(fn,'intensity2','-ascii');
cd('..');

return;
end

function save_tr(tr,x_pos,y_pos,i)

tr;
cd('traces');
fn=['trace_' num2str(x_pos) '-' num2str(y_pos) '-' num2str(i) '.txt'];
save(fn,'tr','-ascii');
cd('..');

return;
end

function combine_molecules

delete('molecules.txt');
fid=fopen('molecules.txt', 'a');
cd('molecules');
files=dir;
numberfiles=length(files);
i=2;
while i < numberfiles,
    i=i+1;
    fid2=fopen(files(i).name,'r');
    if fid2>0
        a=zeros;
        a=dlmread(files(i).name);
        fprintf(fid, '%f\t%f\t%f\t%f\n', a);
        fclose(fid2);
    end
end
fclose(fid);
cd('..');

return;
end

function combine_intensities

delete('intensities.txt');
fid=fopen('intensities.txt', 'a');
cd('intensities');
files=dir;
numberfiles=length(files);
i=2;
while i < numberfiles,
    i=i+1;
    fid2=fopen(files(i).name,'r');
    if fid2>0
        a=zeros;
        a=dlmread(files(i).name);
        fprintf(fid, '%f\n', a);
        fclose(fid2);
    end
end
fclose(fid);
cd('..');

return;
end

function combine_traces

delete('traces.txt');
fid=fopen('traces.txt', 'a');
cd('traces');
files=dir;
numberfiles=length(files);
i=2;
while i < numberfiles,
    i=i+1;
    fid2=fopen(files(i).name,'r');
    if fid2>0
        a=zeros;
        a=dlmread(files(i).name);
        fprintf(fid, '%f\n', a);
        fclose(fid2);
    end
end
fclose(fid);
cd('..');

return;
end