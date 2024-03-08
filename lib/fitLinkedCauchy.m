function [x] = fitLinkedCauchy(horz_x, horzpolar, vert_x, vertpolar)

horznan = ~isnan(horzpolar);
horzpolar = horzpolar(horznan );
horz_x = horz_x(horznan);

vertnan = ~isnan(vertpolar);
vertpolar = vertpolar(vertnan );
vert_x = vert_x(vertnan);


x0 = [1.25 0.25 0.2 0.3 0   1.25 0.25 0.2 0.3 0];
vlb = [0 0 0 0 -1  0 0 0 0 -1];
vub = [3 2 1 1 .3   3 2 1 1 0.3];
options = optimset('fmincon');

x = fmincon(@(x)FitCauchyFunction(x, horz_x, horzpolar, vert_x,vertpolar),x0,[],[],[],[],vlb,vub,@constraintFunc,options);

figure(42); 
clf;
plot(horz_x, horzpolar, vert_x, vertpolar);hold on;
plot(0:0.01:horz_x(end), CauchyLike(x(1:5),0:0.01:horz_x(end)), 'b');
plot(0:0.01:vert_x(end), CauchyLike(x(6:10),0:0.01:vert_x(end)), 'r');
hold off;
drawnow;
    

end



% f = FitCauchyFunction(x, mm_position, horzpolar, vertpolar)
%
function f = FitCauchyFunction(x, horz_x, horzpolar, vert_x, vertpolar)

    horzest = CauchyLike(x(1:5), horz_x);
    vertest = CauchyLike(x(6:10), vert_x);

    
    % Compute fit error as RMSE
    theDiff2 = sum((horzest-horzpolar).^2, 'omitnan') + sum((vertest-vertpolar).^2, 'omitnan');
    f = sqrt( theDiff2/(length(horzpolar)+length(vertpolar)) );

    % figure(42); 
    % clf;
    % plot(horz_x, horzpolar, vert_x, vertpolar);hold on;
    % plot(0:0.01:horz_x(end), CauchyLike(x(1:5),0:0.01:horz_x(end)), 'r');
    % plot(0:0.01:vert_x(end), CauchyLike(x(6:10),0:0.01:vert_x(end)), 'b');
    % hold off;
    % drawnow;
    % pause(0.01);
end

function vals = CauchyLike(x, mm_position)
    vals = x(1)./((1+(mm_position./x(3)).^2)) + x(2)./((1+(mm_position./x(4)).^2)) + x(5);

    % figure(32);    
    % plot(mm_position, x(1)./((1+(mm_position./x(3)).^2)),'b');
    % hold on;
    % plot(mm_position, x(2)./((1+(mm_position./x(4)).^2)),'r');
    % plot(mm_position, vals, 'k');
    % hold off;
    % drawnow;
    % pause(0.01);
end

function [c, ceq] = constraintFunc(x)
    c = [];
    ceq(1) = CauchyLike(x(1:5), 0) - CauchyLike(x(6:10), 0); % Enforce the value at 0 being exactly the same
    ceq(2) = sqrt(sum( (x(1:2) - x(6:7)).^2 )/2); %Enforce the same scale values
end