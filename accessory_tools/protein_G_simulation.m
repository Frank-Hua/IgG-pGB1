%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

This is for simulating trajectories. It assumes both binding and unbinding 
events follow a Poisson process.

Check and adjust parameters that are marked with "frank".
%}

function protein_G_simulation

%initial input
path=command_input('input directory','C:\\Users\\frank\\Documents\\MATLAB','s');
cd(path);

ron=0.008; %frank
roff=0.60; %frank
tr_len=24000; %frank
conc=10; %frank
exposure_time=0.1; %frank
num_mol=272; %frank
num_trans=800; %frank

interarr1=ceil(-10*log(rand(num_trans,num_mol))./ron); %frank
interarr2=ceil(-10*log(rand(num_trans,num_mol))./roff); %frank
total=zeros(3,num_mol);
for i=1:num_mol
    t=1;
    len=sum(interarr1(:,i))+sum(interarr2(:,i));
    tr=zeros(len,1);
    for j=1:num_trans
        tr(t+interarr1(j,i):t+interarr1(j,i)+interarr2(j,i)-1)=1;
        t=t+interarr1(j,i)+interarr2(j,i);
    end
    tr2=tr(1:tr_len); %frank
    total(1,i)=tr_len*conc/nnz(tr2);
    d_tr=conv(tr2,[1,-1],'same');
    d_tr(tr_len)=0;
    index=find(d_tr==1);
    total(2,i)=size(index,1)/(tr_len*exposure_time*conc); %frank
    total(3,i)=total(1,i)*total(2,i); %frank
end

total=total';
fn=[path '\protein G simulation.txt'];
delete(fn);
save(fn,'total','-ascii');

end