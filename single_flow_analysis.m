%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Output files are "molecule_", "intensity_" and "trace_".

    "molecule_" files record background intensity, KM, kon and koff.
    "intensity_" files record intensities when protein G is bound.
    "trace_" files record digital trace information.

For some reason, the pixel size for diffration-limited image is set as 180.0 nm.
The "pixel size" for STORM image is 20.0 nm.

This version is for fast analysis of flow experiment data by
    skipping manual binding site picking.

Check and adjust parameters that are marked with "frank".
%}

function single_flow_analysis

%initial input
%command_input function doesn't do much, but is good to have.
path=command_input('input directory:','C:\\Users\\frank\\Documents\\MATLAB','s');
cd(path);
fname = command_input('input file index #:','1','s');


%read smm files
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
baseline = ones(len,1);
situ = input('movie option: [(0)-256, (1)-flipped 256, (2)-512, (3)-flipped 512] ');

frame = zeros(film_x,film_y,len,'uint16');
for t=1:len
    frame(:,:,t) = fread(fid,[film_x,film_y],'uint16');
end
fclose(fid);


%perform drift correction
%{
"_drift_corrected.txt" files can be directly generated from IDL STORM codes.
    In such a case, use "_drift_corrected.txt" as "_raw.txt" and feed
    drift_correction_control function a mock "_mark.txt" file.
%}
answer = input('skip performing drift correction: [yes-(enter) or no-(no)] ','s');
if strcmp(answer,'no')
    m = dlmread(['film' fname '_raw.txt']);
    m(:,5) = uint16(m(:,5));
    m0 = dlmread(['film' fname '_mark.txt']);
    m0(:,5) = uint16(m0(:,5));
    drift_correction_control(fname,m,m0);
end

%{
"_drift_corrected_histogram.tif" is basically a STORM image of all binding sites.
    Its pixel number is 81=9*9 times of the diffration-limited image.
"_drift.txt" is a drift file. Drifts in x and y are recorded in the unit of nm.
%}
m2 = dlmread(['film' fname '_drift_corrected.txt']);
m2(:,5) = uint16(m2(:,5));
m3 = imread(['film' fname '_drift_corrected_histogram.tif'],'TIFF');
m4 = dlmread(['film' fname '_drift.txt']);


%plot spot number per frame
%This can be used to calculate photobleaching rates.
answer2 = input('skip calculating spot number per frame: [yes-(enter) or no-(no)] ','s');
if strcmp(answer2,'no')
    spot_number_per_frame(fname,m2);    
end


n = 0;

%correct trace baseline
answer3 = input('skip correcting trace baseline: [yes-(enter) or no-(no)] ','s');
if strcmp(answer3,'no')
    hdl = figure;
    temp = dlmread(['film' fname '_baseline_positions.txt']);
    center_xb=temp(:,1);
    center_yb=temp(:,2);
    for i=1:length(center_xb)
        t_baseline = baseline_localization(n,m2,center_xb(i),center_yb(i),frame,len,m4,hdl,situ);
        baseline = [baseline t_baseline];
        if size(baseline,2) > 1
            baseline(:,1) = mean(baseline(:,2:end),2);
        end
    end
end


answer4 = input('specify binding site finding parameter-circle radius in sr-pixels: [default=5] ');
rad=num2str(answer4);

%Boundaries are from 1 to 2295 super-resolution (sr)-pixels.
answer5 = command_input('specify binding site finding parameter-x lower boundary in sr-pixels:','206','');
answer6 = command_input('specify binding site finding parameter-x upper boundary in sr-pixels:','2090','');
answer7 = command_input('specify binding site finding parameter-y lower boundary in sr-pixels:','206','');
answer8 = command_input('specify binding site finding parameter-y upper boundary in sr-pixels:','2090','');

x1=num2str(answer5);
x2=num2str(answer6);
y1=num2str(answer7);
y2=num2str(answer8);

%{
Analyzed binding site will be deleted from temp_good. That way we keep a
    record of what has been analyzed.
%}
fn = ['film' fname '_localizations_' x1 '_' x2 '_' y1 '_' y2 '_' rad '.txt'];
if exist(fn)
    good=dlmread(fn);
    no_good=size(good,1);
else
    if isempty(answer4)
        [good,no_good] = finding_localization(m3,answer5,answer6,answer7,answer8,2000);
    else
        [good,no_good] = finding_localization_radius(m3,answer5,answer6,answer7,answer8,2000,answer4);
    end
    [good,no_good] = redundant_binding_site(path,good,no_good);
end
disp(['total number ' num2str(no_good)]);


hdl2 = figure;
hdl3 = figure;

mkdir('molecules');
mkdir('intensities');
mkdir('traces');
index = [];

while n < no_good
    n=n+1;
    %Convert good(n,:) into center_xy with values in the unit of nm.
    center_x = good(n,1)*20.0;
    center_y = good(n,2)*20.0;
    [total,intensity2,tr,value]=analyze_localization(n,m2,center_x,center_y,frame,len,m4,baseline,hdl2,hdl3,situ);
    if value == -1
        save_molecule(total,floor(good(n,1)),floor(good(n,2)),n);
        save_intensity2(intensity2,floor(good(n,1)),floor(good(n,2)),n);
        save_tr(tr,floor(good(n,1)),floor(good(n,2)),n);
        index = [index n];
    end
    if value == -2
        break;
    end
    if value > 0
        n=value-1;
    end
end

if ~isempty(index)
    good(1:index(end),:) = [];
end
close('all');
combine_molecules;
combine_intensities
combine_traces;
save_localization(fn,good);

end







function save_localization(fn,good)

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