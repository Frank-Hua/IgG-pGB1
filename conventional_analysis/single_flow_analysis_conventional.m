%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Output files are "molecule_xx-yy-#", "intensity_xx-yy-#" and "trace_xx-yy-#".

    "molecule_xx-yy-#" files record background intensity, KM, kon and koff.
    "intensity_xx-yy-#" files record intensities when protein G is bound.
    "trace_xx-yy-#" files record the digital trace.

For some reason, the pixel size for diffration-limited image is set as 180.0 nm. And the "pixel size" for STORM image is 20.0 nm.

This version is for fast analysis by skipping manual binding site picking.

This version is for conventional analysis.

Adjustable parameters are marked with "frank".
%}

function single_flow_analysis_conventional

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


%% correct trace baseline
answer3 = input('skip correcting intensity baseline: [yes-(enter) or no-(no)] ','s');
if strcmp(answer3,'no')
    temp = dlmread(['film' fname '_baseline_positions.txt']);
    center_xb=temp(:,1);
    center_yb=temp(:,2);
    hdl = figure;
    baseline = ones(len,1);
    for i=1:length(center_xb)
        t_baseline = baseline_site(0,m2,center_xb(i),center_yb(i),frame,len,m4,hdl,situ);
        baseline = [baseline t_baseline];
        if size(baseline,2) > 1
            baseline(:,1) = mean(baseline(:,2:end),2);
        end
    end
    delete(['film' fname '_baseline.txt']);
    fid = fopen(['film' fname '_baseline.txt'], 'a' );
    for j=1:size(baseline(:,1),1)
        fprintf(fid,'%f\n',baseline(j,1));
    end
    fclose(fid);
end

baseline = dlmread(['film' fname '_baseline.txt']);


%% start analysis
mkdir('molecules');
mkdir('intensities');

answer4 = command_input('specify binding site finding parameter-circle radius in sr-pixels:','14','');
rad=num2str(answer4);

%Boundaries are between 1 to 2295 super-resolution (sr)-pixels.
answer5 = command_input('specify binding site finding parameter-x lower boundary in sr-pixels:','206','');
answer6 = command_input('specify binding site finding parameter-x upper boundary in sr-pixels:','2090','');
answer7 = command_input('specify binding site finding parameter-y lower boundary in sr-pixels:','206','');
answer8 = command_input('specify binding site finding parameter-y upper boundary in sr-pixels:','2090','');

x1=num2str(answer5);
x2=num2str(answer6);
y1=num2str(answer7);
y2=num2str(answer8);

answer9 = command_input('specify binding site finding parameter-minimum localization event number:','100','');
%For some reason, this value needs to be multiplied by 20.
answer9 = answer9*20.0;
thres=num2str(answer9);

%{
Analyzed binding site will be deleted from good. That way we keep a record of what has been analyzed.
%}
fn = ['film' fname '_sites_' x1 '_' x2 '_' y1 '_' y2 '_' rad '_' thres '.txt'];
if exist(fn,'file')
    good=dlmread(fn);
    no_good=size(good,1);
else
    [good,no_good] = finding_site_radius_conventional(m2,answer5,answer6,answer7,answer8,answer9,answer4);
    [good,no_good] = redundant_binding_sites(path,good,no_good);
end
disp(['total number ' num2str(no_good)]);

n = 0;
hdl2 = figure;
hdl3 = figure;
index = [];

while n < no_good
    n=n+1;
    %Convert good(n,:) into center_xy with values in the unit of nm.
    center_x = good(n,1)*20.0;
    center_y = good(n,2)*20.0;
    [total,intensity2,value]=analyze_site_conventional(n,m2,center_x,center_y,frame,len,m4,baseline,hdl2,hdl3,situ);
    if value == -1
        save_molecule(total,floor(good(n,1)),floor(good(n,2)),n);
        save_intensity2(intensity2,floor(good(n,1)),floor(good(n,2)),n);
        index = [index n];
    end
    if value == -2
        break;
    end
    if value > 0
        n=value-1;
    end
end

index=sort(index,'ascend');
if ~isempty(index)
    good(1:index(end),:) = [];
end

close('all');
save_site(fn,good);
combine_molecules;
combine_intensities

end







function save_site(fn,good)

good;
delete(fn);
save(fn,'good','-ascii');

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