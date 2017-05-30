%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.
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