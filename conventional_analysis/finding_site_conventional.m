%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function [good,no_good]=finding_site_conventional(m,n1,n2,n3,n4,threshold)

r=14;
[good,no_good]=finding_site_radius_conventional(m,n1,n2,n3,n4,threshold,r);

return;
end