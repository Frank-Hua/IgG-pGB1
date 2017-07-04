%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function circle = circle_generator(r)

circle = zeros(2*r+1,2*r+1);
for i = 1:(2*r+1)
    for j = 1:(2*r+1)
        dist = sqrt((i-r-1)^2+(j-r-1)^2);
        if round(dist) == r
            circle(i,j)=1;
        end
    end
end

end