%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Adjustable parameters are marked with "frank".
%}

function value=command_input(title,default,format)

switch format
    case 's'
        value = input([title ' [default=' default '] '],'s');
        if isempty(value)
            value = default;
        end
    otherwise
        value = input([title ' [default=' default '] ']);
        if isempty(value)
            value = str2num(default);
        end
end

return;
end