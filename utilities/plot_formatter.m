%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function temp=plot_formatter(title_s,p1,p2,p3,p4,x1,x2,y1,y2)

if ~isempty(title_s)
    title(title_s);
end

if p1;  grid on;  end
if p2;  zoom on;  end

axis tight;
temp=axis;
if x1 ~= 'off';	temp(1)=x1;	end
if x2 ~= 'off';	temp(2)=x2;	end
if y1 ~= 'off';	temp(3)=y1;	end
if y2 ~= 'off';	temp(4)=y2; end
axis(temp);

if ~p3;  axis off;   end
if p4;  hold on;    end

end