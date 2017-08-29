%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Adjustable parameters are marked with "frank".
%}

function s_avg_dist2=flip_drift_correction(s_avg_dist,flip)

if flip
    s_avg_dist2(:,2)=s_avg_dist(:,1);
    s_avg_dist2(:,1)=s_avg_dist(:,2);
else
    s_avg_dist2=s_avg_dist;
end

end