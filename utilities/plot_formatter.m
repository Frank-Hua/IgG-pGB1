%{
This is made for analyzing the kinetics data of protein G-IgG interaction, a project in collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Adjustable parameters are marked with "frank".
%}

function temp=plot_formatter(title_s,p1,p2,p3,p4,x1,x2,y1,y2)

if ~isempty(title_s)
    title(title_s);
end

if p1;  grid on;  end
if p2;  zoom on;  end

axis tight;
temp=axis;
if ~strcmp(x1,'off');	temp(1)=x1;	end
if ~strcmp(x2,'off');	temp(2)=x2;	end
if ~strcmp(y1,'off');	temp(3)=y1;	end
if ~strcmp(y2,'off');	temp(4)=y2; end
axis(temp);

if ~p3;  axis off;   end
if p4;  hold on;    end

end