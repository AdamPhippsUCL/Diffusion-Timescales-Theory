clear;

R = 10;


f=figure;
f.Position = [680   458   500   500];
ax = axes;
ax.Position = [0.1, 0.1, 0.8, 0.8];

% Plot circle
thetas = 0:0.02:2*pi;
circx = R*cos(thetas);
circy = R*sin(thetas);
plot(circx, circy, color = 'k')
hold on
xline(0, '--', color=[.2 .2 .2]);
yline(0, '--', color=[.2 .2 .2]);

xlim([-12, 12])
ylim([-12, 12])



% Simulate diffusion

stepsize = 0.2;
Nstep = 5000;

positions = zeros(Nstep, 2);

CData = zeros(Nstep, 3);
cols = colororder;
cindx = 1;

startindx = 1;

beta = 2*pi*rand(1);
f=0.8;

for stepindx = 1:Nstep-1

    alpha = 2*pi*rand(1);

    vec = [cos(alpha), sin(alpha)] + f*[cos(beta), sin(beta)];

    step = stepsize*(1/norm(vec))*vec;

    positions(stepindx+1, :) = positions(stepindx,:) + step;

    new_pos = [positions(stepindx+1, :)];


    if new_pos(1)^2 + new_pos(2)^2 >= R^2
        step = -step;
        positions(stepindx+1, :) = positions(stepindx,:) + step;


        if norm(positions(stepindx,:) - positions(startindx,:)) > 8

            cindx = cindx+1;
        
            if cindx==8
                return
            end
    
            plot(positions(startindx:stepindx,1), positions(startindx:stepindx,2), color = [cols(cindx,:)])
    
            startindx = stepindx;

      
        end

        beta = 2*pi*rand(1);

    end


    if randi(50)==1
        beta = 2*pi*rand(1);
    end


    % CData(stepindx,:)=cols(cindx,:);

end


% plot(positions(:,1), positions(:,2), '--',color=[.2 .2 .2]);
% scatter(positions(:,1), positions(:,2), '.',CData=CData);