%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Adjustable parameters are marked with "frank".
%}

function circle = circle_generator(r)

circle = zeros(2*r+1,2*r+1);
for i = 1:(2*r+1)
    for j = 1:(2*r+1)
        dist = sqrt((i-r-1)^2+(j-r-1)^2);
        if round(dist) == r-1
            circle(i,j)=1;
        end
    end
end

end