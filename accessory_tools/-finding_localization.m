%{
This is for analyzing the protein G-IgG binding kinetics, a project in
    collaboration with Prof. Wei Cheng in UMich, Ann Arbor.

Check and adjust parameters that are marked with "frank".
%}

function [good,no_good]=finding_localization(count,n1,n2,n3,n4,threshold)

circle=zeros(11,11);
circle(:,1) = [ 0,0,0,0,0,0,0,0,0,0,0];
circle(:,2) = [ 0,0,0,0,1,1,1,0,0,0,0];
circle(:,3) = [ 0,0,0,1,0,0,0,1,0,0,0];
circle(:,4) = [ 0,0,1,0,0,0,0,0,1,0,0];
circle(:,5) = [ 0,1,0,0,0,0,0,0,0,1,0];
circle(:,6) = [ 0,1,0,0,0,0,0,0,0,1,0];
circle(:,7) = [ 0,1,0,0,0,0,0,0,0,1,0];
circle(:,8) = [ 0,0,1,0,0,0,0,0,1,0,0];
circle(:,9) = [ 0,0,0,1,0,0,0,1,0,0,0];
circle(:,10) = [ 0,0,0,0,1,1,1,0,0,0,0];
circle(:,11) = [ 0,0,0,0,0,0,0,0,0,0,0];
g_peaks = zeros(3,3,7,7);
for k = 1:3
    for l = 1:3
        offx = 0.5*(k-2);
        offy = 0.5*(l-2);
        for i = 1:7
            for j = 1:7
                dist = 0.69*((i-4+offx)^2+(j-4+offy)^2); %frank
                g_peaks(k,l,i,j) = exp(-dist);
            end
        end
    end
end

good = zeros(8000,2);
foob = zeros(49,1);
diff = zeros(3,3);
no_good = 0;
for i = n1:n2
    for j = n3:n4
        if count(i,j) > 0
            dum = count(i-3:i+3,j-3:j+3);
            foob = dum(:);
            sum_foob = sum(foob);
            [z,foo] = max(foob);
            x = mod((foo-1),7)-3;
            y = floor((foo-1)/7)-3;
            if (x == 0) && (y == 0)
                x = x+i;
                y = y+j;
                quality = 1;
                for k = -5:5
                    for l = -5:5
                        if circle(k+6,l+6) > 0
                            if count(x+k,y+l) > 0.45*z || sum_foob <= threshold %thresholds to select localizations
                                quality = 0;
                            end
                        end
                    end
                end

                if quality == 1
                    cur_best = 1000000.0;
                    for k = 1:3
                        for l = 1:3
                            dum2 = reshape(g_peaks(k,l,:,:),7,7);
                            dum3 = sum(sum(abs(double(z)*dum2(:,:)-double(count(x-3:x+3,y-3:y+3)))));
                            diff(k,l) = dum3;
                            if diff(k,l) < cur_best
                                best_x = k;
                                best_y = l;
                                cur_best = diff(k,l);
                            end
                        end
                    end
                    
                    fit_x = double(x)-0.5*(best_x-2);
                    fit_y = double(y)-0.5*(best_y-2);
                    no_good = no_good+1;
                    good(no_good,1) = fit_x;
                    good(no_good,2) = fit_y;
                end
            end
        end
    end
end

good = good(1:no_good,:);

return;
end