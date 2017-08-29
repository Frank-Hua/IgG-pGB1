%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Adjustable parameters are marked with "frank".
%}

function [c_intensity,baseline2] = baseline_correction(intensity,tr,bsl0,bsl)

c_intensity = intensity;

len = length(tr);
if len < 5
    warning('total length of this movie is too short');
end

%Use value 10 to expand protein G binding peaks wider.
tr2 = ones(len,1)-imdilate(tr,strel('line',10,90)); %frank
intensity2 = intensity.*tr2;

baseline(1) = bsl0;
for i = 2 : floor(len/5)
    temp = nonzeros(intensity2((i-1)*5+1:i*5));
    if ~isempty(temp)
        baseline(i) = mean(temp);
    else
        baseline(i) = baseline(i-1);
    end
end
baseline(end+1) = baseline(end);

x = 1:5:floor(len/5)*5+1;
xx = 1:floor(len/5)*5;
baseline2 = pchip(x,baseline,xx)';

% baseline2 = medfilt1(baseline2,51);
baseline2 = smooth(baseline2,51);

c_intensity(1:floor(len/5)*5) = intensity(1:floor(len/5)*5)./baseline2*bsl;
    
end