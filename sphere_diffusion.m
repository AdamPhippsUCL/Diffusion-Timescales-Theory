clear;

%% Settings

gamma=2.675; % *1e8 s^-1 T^-1;
D=2;
R=10;

tau = R^2/(2*D);

% Lobe area
A=200;

% b-value
b=1000;


%% Change with DELTA

f=figure;
f.Position  = [286.6 97 916.8 657.6];
ax = axes;

h = imagesc(ax, [0 2*tau], [0 1], linspace(0,2*tau,256)');
cmap=colormap(ax, 'turbo');
caxis([0 1])

cb = colorbar(ax);
cb.Ticks = [0, 0.5, 1];
cb.TickLabels = {'0', '$\Delta=\tau$', '$\Delta=2\tau$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 18;
set(h, 'Visible','off')
axis xy
hold on



deltamin = 0.1;
deltamax = 51;
deltastep = 1;
deltas = deltamin:deltastep:deltamax;

Deltamin = deltamin;
Deltastep = 1;
Deltamax = 80;
Deltas = Deltamin:Deltastep:Deltamax;
NDelta = length(Deltas);

phase_vars = zeros(length(deltas), length(Deltas));
phase_vars1 = zeros(length(deltas), 1);

for dindx = 1:length(deltas)

    delta = deltas(dindx);
    thisDeltas = [delta:Deltastep:Deltamax];
    N=length(thisDeltas);


    % Fixed lobe area
    G=A/delta;
    
    for indx = 1:N

        % G = stejskal(delta, thisDeltas(indx), bval=1000);
        sphere_signal = sphereGPD(delta, thisDeltas(indx), G, R, D);
        
        % b= stejskal(delta, thisDeltas(indx), G=G);
        % sphere_signal = ball(b, D);

        phase_vars(dindx, (NDelta-N)+indx) = -2*log(sphere_signal);

        if indx == 1
            phase_vars1(dindx) = phase_vars(dindx, (NDelta-N)+indx);
        end
    
    end

    plot(thisDeltas, phase_vars(dindx, (NDelta-N)+1:NDelta), Color=cmap(ceil(256*delta/deltamax), :), LineWidth = 1, HandleVisibility='off')
    hold on

    if dindx == 1
            
        % alpha = phase_vars(dindx,1)/(Deltas(1)-delta/3);
        alpha = 2*((gamma*1e8)^2)*((A*1e-6)^2)*(2*1e-9)*(1e-3);  % 2*gamma^2*A^2*D
               
        plot(Deltas, alpha*(thisDeltas-delta/3), '--', color = 'k', LineWidth = 1, DisplayName='$\Delta - \frac{\delta}{3}$', HandleVisibility='off')
    
        % Compute saturation value
        phase_var_sat = 2*((gamma*1e8)^2)*((A*1e-6)^2)*( 0.2*((R*1e-6)^2) ); %- (D*1e-9)*(delta*1e-3)/3)
        yline(phase_var_sat, '--', color = 'k', LineWidth = 1, HandleVisibility='off')
    end

end

plot(deltas, phase_vars1(:), '--', color = 'k', HandleVisibility='off')

[~, indx] = min(abs(deltas-tau));
plot([tau, tau], [0, phase_vars1(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')

[~, indx] = min(abs(deltas-2*tau));
plot([2*tau, 2*tau], [0, phase_vars1(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')

c=colorbar;
c.Ticks = [0 0.5 1];
c.TickLabels = {'$\delta=0$', '$\delta = \tau$', '$\delta = 2\tau$'} ;
c.set('FontSize', 18)
c.TickLabelInterpreter = 'latex';

xlim([0, max(Deltas)])
ylim([0, 0.15])

yticks([0, phase_var_sat])
yticklabels({'', '$\frac{2}{5}\gamma^2 A^2 R^2$'})
xticks([0, tau, 2*tau])
xticklabels({'0', '$\tau$', '$2\tau$'})
ax = gca();
set(ax, 'TickLabelInterpreter', 'latex')
ax.FontSize = 20;
ax.YAxis.FontSize = 18;

text(14.0, 0.135, '$2 \gamma^2 A^2 D ( \Delta - \frac{\delta}{3})$', Interpreter='latex', FontSize=18)

xlabel('$\Delta \rightarrow$', Position=[max(Deltas)-4,-0.002],Interpreter='latex', FontSize=20)
ylabel('$\langle \Delta\phi^2 \rangle \rightarrow$', Position=[-1 0.144] , Interpreter='latex', FontSize=20)

saveas(f, 'sphere_phase_var_big_DELTA.png')



figure;

% Specific scan
delta = 10;
Delta = 35;
G=A/delta;
bval = stejskal(delta, Delta, G=G);

sphere_signal = sphereGPD(delta, Delta, G, R, D);
pv_sphere = -2*log(sphere_signal);

scatter(Delta, exp(-pv_sphere/2))
hold on

ball_signal = ball(bval, D);
pv_ball = -2*log(ball_signal);
scatter(Delta, exp(-pv_ball/2))


frac = pv_sphere/pv_ball;

%% Look at change with delta (plot <phi^2> + 2/3 gamma^2 A^2 D delta)

f1=figure;
f1.Position  = [286.6 97 916.8 657.6];
ax1 = axes;

this_DELTA_max = 56;

h = imagesc(ax1, [0 2*tau], [0 1], linspace(0,2*tau,256)');
colormap(ax1, 'turbo')
caxis([0 this_DELTA_max])

cb = colorbar(ax1);
cb.Ticks = [0, tau, 2*tau];
cb.TickLabels = {'0', '$\Delta=\tau$', '$\Delta=2\tau$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 18;
set(h, 'Visible','off')
axis xy
hold on

inds = [1:1:this_DELTA_max];
last_points = [];

% Predicted gradient
pred_grad = (-2/3)*((gamma*1e8)^2)*((A*1e-6)^2)*(2*1e-9)*(1e-3);

for Delta =Deltas(inds)

    [~,Deltaindx] = min(abs(Deltas-Delta));
    
    this_phase_vars = phase_vars(:,Deltaindx);
    bool = this_phase_vars~=0;
    % beta = phase_vars(1,1)/(Deltamin-deltamin/3);
    % beta=alpha;

    this_deltas = deltas(bool);
    this_y = this_phase_vars(bool)-pred_grad*deltas(bool)';

    % Plot
    plot(ax1, [0 this_deltas], [this_y(1); this_y], Color = cmap(min([ceil(256*Delta/(this_DELTA_max)), 256]), :), LineWidth = 1);
    hold(ax1, 'on');

    last_points = [last_points, this_y(end)];

    
end

yline(phase_var_sat+pred_grad*0.1, '--', LineWidth = 1, color = 'k')
plot([0, Deltas(inds)], [0, last_points], '--', color = 'k')

[~, indx] = min(abs(Deltas(inds)-tau));
plot([tau, tau], [0, last_points(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')


% plot([0 Delta], [last_points(end)+pred_grad*2*tau last_points(end)])

xlim([0 50])
ylim([0 0.25])

yticks([0, phase_var_sat])
yticklabels({'', '$\frac{2}{5}\gamma^2 A^2 R^2$'})
xticks([0, tau])
xticklabels({'0', '$\tau$'})
set(ax1, 'TickLabelInterpreter', 'latex')
ax1.FontSize = 20;
ax1.YAxis.FontSize = 18;

xlabel('$\delta \rightarrow$', Position=[48,-0.002],Interpreter='latex', FontSize=20)
ylabel('$\langle \Delta\phi^2 \rangle + \frac{2}{3} \gamma^2 A^2 D \delta \rightarrow$', Position=[-0.64 0.204] , Interpreter='latex', FontSize=18)
saveas(f1, 'sphere_phase_var_small_delta_1.png')



%% Look at change with delta (plot <phi^2>)

f2=figure;
f2.Position  = [286.6 97 916.8 657.6];
ax2 = axes;

this_DELTA_max = 56;

h = imagesc(ax2, [0 2*tau], [0 1], linspace(0,2*tau,256)');
colormap(ax2, 'turbo')
caxis([0 this_DELTA_max])

cb = colorbar(ax2);
cb.Ticks = [0, tau, 2*tau];
cb.TickLabels = {'0', '$\Delta=\tau$', '$\Delta=2\tau$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 18;
set(h, 'Visible','off')
axis xy
hold on

inds = [1:1:this_DELTA_max];
last_points = [];

% Predicted gradient
pred_grad = (-2/3)*((gamma*1e8)^2)*((A*1e-6)^2)*(2*1e-9)*(1e-3);

for Delta =Deltas(inds)

    [~,Deltaindx] = min(abs(Deltas-Delta));
    
    this_phase_vars = phase_vars(:,Deltaindx);
    bool = this_phase_vars~=0;
    % beta = phase_vars(1,1)/(Deltamin-deltamin/3);
    % beta=alpha;

    this_deltas = deltas(bool);
    this_y = this_phase_vars(bool);

    % Plot
    plot(ax2, [0 this_deltas], [this_y(1); this_y], Color = cmap(min([ceil(256*Delta/(this_DELTA_max)), 256]), :), LineWidth = 1);
    hold(ax2, 'on');

    last_points = [last_points, this_y(end)];

    
end

yline(phase_var_sat+pred_grad*0.1, '--', LineWidth = 1, color = 'k')
plot([0, Deltas(inds)], [0, last_points], '--', color = 'k')

[~, indx] = min(abs(Deltas(inds)-tau));
plot([tau, tau], [0, last_points(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')


% plot([0 Delta], [last_points(end)+pred_grad*2*tau last_points(end)])

xlim([0 50])
ylim([0 0.15])

yticks([0, phase_var_sat])
yticklabels({'', '$\frac{2}{5}\gamma^2 A^2 R^2$'})
xticks([0, tau])
xticklabels({'0', '$\tau$'})
set(ax2, 'TickLabelInterpreter', 'latex')
ax2.FontSize = 20;
ax2.YAxis.FontSize = 18;

xlabel('$\delta \rightarrow$', Position=[48,-0.002],Interpreter='latex', FontSize=20)
ylabel('$\langle \Delta\phi^2 \rangle \rightarrow$', Position=[-0.64 0.144] , Interpreter='latex', FontSize=18)
saveas(f2, 'sphere_phase_var_small_delta_2.png')