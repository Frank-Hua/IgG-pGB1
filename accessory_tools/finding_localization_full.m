%This is for analyzing the protein G-IgG binding kinetics, a project in
%collaboration with Prof. Wei Cheng in UMich, Ann Arbor.
%The "_full" version is for finding localizations in the whole area.
%Changable parameters are labeled with "frank".

function finding_localization_full

path=command_input('input directory','C:\\Users\\frank\\Documents\\MATLAB','s');
cd(path);
fname = command_input('input file index #','1','s');
count = imread(['film' fname '_drift_corrected_histogram.tif'],'TIFF');

[good,no_good]=finding_localization(count,206,2090,206,2090,2000); %frank
%[good,no_good]=finding_localization(count,673,1967,495,1779,2000);

fn=['film' fname '_localizations_full.txt'];
save_localization(fn,good);
disp(no_good);

end







function save_localization(fn,good)

good;
delete(fn);
save(fn,'good','-ascii');

return;
end