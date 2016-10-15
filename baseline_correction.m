%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function [c_intensity,baseline2] = baseline_correction(intensity,tr,bsl0,bsl)

c_intensity = intensity;

len = length(tr);

%10 expands protein G binding peaks wider.
tr2 = ones(len,1)-imdilate(tr,strel('line',10,90)); %frank
intensity2 = intensity.*tr2;

baseline = ones(floor(len/5)+1,1);
baseline(1) = bsl0;
for i = 2 : floor(len/5)
    temp = nonzeros(intensity2((i-1)*5+1:i*5));
    if ~isempty(temp)
        baseline(i) = min(temp);
    else
        baseline(i) = baseline(i-1);
    end
end
baseline(end) = baseline(end-1);

x = 1:5:floor(len/5)*5+1;
xx = 1:floor(len/5)*5;
baseline2 = pchip(x,baseline,xx)';

c_intensity(1:floor(len/5)*5) = intensity(1:floor(len/5)*5)./baseline2*bsl;
    
end