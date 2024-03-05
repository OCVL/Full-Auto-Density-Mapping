function [x] = fitLinkedCauchy(horz_x, horzpolar, vert_x, vertpolar)

horznan = ~isnan(horzpolar);
horzpolar = horzpolar(horznan );
horz_x = horz_x(horznan);

vertnan = ~isnan(vertpolar);
vertpolar = vertpolar(vertnan );
vert_x = vert_x(vertnan);


x0 = [1.75 0.4 0.2 0.3  1.75 0.4 0.2 0.3];
vlb = [];
vub = [];
options = optimset('fmincon');

x = fmincon(@(x)FitCauchyFunction(x, horz_x, horzpolar, vert_x,vertpolar),x0,[],[],[],[],vlb,vub,@constraintFunc,options);

figure(42); 
clf;
plot(horz_x, horzpolar, vert_x, vertpolar);hold on;
plot(0:0.01:horz_x(end), CauchyLike(x(1:4),0:0.01:horz_x(end)), 'b');
plot(0:0.01:vert_x(end), CauchyLike(x(5:8),0:0.01:vert_x(end)), 'r');
hold off;
drawnow;
    

end



% f = FitCauchyFunction(x, mm_position, horzpolar, vertpolar)
%
function f = FitCauchyFunction(x, horz_x, horzpolar, vert_x, vertpolar)

    horzest = CauchyLike(x(1:4), horz_x);
    vertest = CauchyLike(x(5:8), vert_x);

    
    % Compute fit error as RMSE
    theDiff2 = sum((horzest-horzpolar).^2, 'omitnan') + sum((vertest-vertpolar).^2, 'omitnan');
    f = sqrt( theDiff2/(length(horzpolar)+length(vertpolar)) );

    % figure(42); 
    % clf;
    % plot(horz_x, horzpolar, vert_x, vertpolar);hold on;
    % plot(0:0.01:horz_x(end), CauchyLike(x(1:4),0:0.01:horz_x(end)), 'r');
    % plot(0:0.01:vert_x(end), CauchyLike(x(5:8),0:0.01:vert_x(end)), 'b');
    % hold off;
    % drawnow;
    % pause(0.01);
end

function vals = CauchyLike(x, mm_position)
    vals = x(1)./((1+(mm_position./x(3)).^2)) + x(2)./((1+(mm_position./x(4)).^2));
end

function [c, ceq] = constraintFunc(x)
    c = [];
    ceq(1) = CauchyLike(x(1:4), 0) - CauchyLike(x(5:8), 0); % Enforce the value at 0 being exactly the same
    ceq(2) = sqrt(sum( (x(1:2) - x(5:6)).^2 )/2); %Enforce the same scale values
end