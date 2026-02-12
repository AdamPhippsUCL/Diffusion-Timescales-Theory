% New sphere diffusion simulation script

% Separate simulation over grid and figure creation into different
% sections.


clear;


%% Settings

gamma=2.675; % *1e8 s^-1 T^-1;

% Diffusivity
D=2; % *1e-9 m^2 s^-1

% Radius
R=10; % *1e-6 m

% Gradient lobe area
A = 200; % *1 mT ms m^-1 = *1e-6 Ts m^-1

% Tau timescale
tau = R^2/(2*D);

% Phase variance saturation (prediction from theory)
phase_var_sat = 2*((gamma*1e8)^2)*((A*1e-6)^2)*( 0.2*((R*1e-6)^2) ); % 2 gamma^2 A^2 R^2/5

% Phase variance linear gradient (prediction from theory)
alpha = 2*((gamma*1e8)^2)*((A*1e-6)^2)*(2*1e-9)*(1e-3); % 2*gamma^2*A^2*D  


% == Grid

% delta
deltamin = 0.1;
deltamax = 501;
deltastep = 1;
deltas = deltamin:deltastep:deltamax;
Ndelta = length(deltas);

% Delta
Deltamin = deltamin;
Deltastep = 1;
Deltamax = 510;
Deltas = Deltamin:Deltastep:Deltamax;
NDelta = length(Deltas);

% Initialise phase variance array
phase_vars = zeros(Ndelta, NDelta);





%% Simulate over grid

for dindx = 1:Ndelta

    delta = deltas(dindx);
    % thisDeltas = Deltas(Deltas>delta);
    thisDeltas = [delta:Deltastep:Deltamax];
    thisNDelta = length(thisDeltas);

    % Gradient strength (fixed lobe area)
    thisG = A/delta;

    for Dindx = 1:thisNDelta
        
        % Simulate sphere phase variance
        sphere_signal = sphereGPD(delta, thisDeltas(Dindx), thisG, R, D);
        phase_vars(dindx, (NDelta-thisNDelta)+Dindx) = -2*log(sphere_signal);

    end

end



% Showing 1/delta dependence
this_ys = phase_vars(:,end).*deltas';
f=figure;
f.Position  = [286.6 97 916.8 657.6];
ax = axes;
ax.FontSize = 20;
plot(deltas, this_ys, LineWidth=1.2, DisplayName=['\Delta = ' num2str(floor(Deltamax/tau)) '\tau'])
xlim([0 450])
xticks([0,4*tau,8*tau,  12*tau])
xticklabels({ '0', '4\tau',  '8\tau',  '12\tau'})
ylim([0 3.5])
yticks([0, 1.005*max(this_ys)])
yticklabels({'', ''})
ax.FontSize = 18;
xlabel('$\delta \rightarrow$', Position=[430,-0.04],Interpreter='latex', FontSize=22)
ylabel('$\langle \Delta\phi^2 \rangle \cdot \delta \rightarrow$', Position=[-7 3.2] , Interpreter='latex', FontSize=20)
yline(1.005*max(this_ys), '--', LineWidth=1, color = .2*[1,1,1], HandleVisibility='off')
legend;
saveas(f, 'sphere_phase_var_1_over_delta_dependence.png')




%% Plot phase variance as function of DELTA

% Initialise figure
f1=figure;
f1.Position  = [286.6 97 916.8 657.6];
ax1 = axes;

% Create color bar
h = imagesc(ax1, [0 2*tau], [0 1], linspace(0,2*tau,256)');
cmap=colormap(ax1, 'turbo');
caxis([0 1])
cb = colorbar(ax1);
cb.Ticks = [0, 0.5, 1];
cb.TickLabels = {'0', '$\delta=\tau$', '$\delta=2\tau$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 18;
set(h, 'Visible','off')
axis xy
hold on

phase_vars_firsts = zeros(Ndelta, 1);

for dindx = 1:Ndelta

    delta = deltas(dindx);
    thisDeltas = [delta:Deltastep:Deltamax];
    thisNDelta = length(thisDeltas);

    plot(thisDeltas, phase_vars(dindx, (NDelta-thisNDelta)+1:NDelta), Color=cmap(ceil(256*delta/deltamax), :), LineWidth = 1, HandleVisibility='off')

    phase_vars_firsts(dindx) = phase_vars(dindx, (NDelta-thisNDelta)+1);; 

end

% Connect start points of each curve
plot(deltas, phase_vars_firsts, '--', color = 'k', HandleVisibility='off');

% Mark tau and 2*tau
[~, indx] = min(abs(deltas-tau));
plot([tau, tau], [0, phase_vars_firsts(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')

[~, indx] = min(abs(deltas-2*tau));
plot([2*tau, 2*tau], [0, phase_vars_firsts(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')


% Phase variance saturation
yline(phase_var_sat, '--', color = 'k', LineWidth = 1, HandleVisibility='off')

% Phase variance linear growth (delta, Delta << tau)
plot(Deltas, alpha*(Deltas-deltas(1)/3), '--', color = 'k', LineWidth = 1, HandleVisibility='off')
text(14.0, 0.135, '$2 \gamma^2 A^2 D ( \Delta - \frac{\delta}{3})$', Interpreter='latex', FontSize=18)

% Axis limits
xlim([0, max(Deltas)])
ylim([0, 0.15])

% Axis ticks
yticks([0, phase_var_sat])
yticklabels({'', '$\frac{2}{5}\gamma^2 A^2 R^2$'})
xticks([0, tau, 2*tau])
xticklabels({'0', '$\tau$', '$2\tau$'})
set(ax1, 'TickLabelInterpreter', 'latex')
ax1.FontSize = 20;
ax1.YAxis.FontSize = 18;

% Axis labels
xlabel('$\Delta \rightarrow$', Position=[max(Deltas)-3,-0.002],Interpreter='latex', FontSize=20)
ylabel('$\langle \Delta\phi^2 \rangle \rightarrow$', Position=[-1 0.144] , Interpreter='latex', FontSize=20)


saveas(f1, 'sphere_phase_var_big_DELTA.png')


%% Plot as function of delta

thisDeltamax = 56; % Maximum DELTA value to plot for
Imax = find(Deltas>thisDeltamax, 1, 'first');
inds = 1:Imax;


f2=figure;
f2.Position  = [286.6 97 916.8 657.6];
ax2 = axes;

% Create colorbar
h = imagesc(ax2, [0 2*tau], [0 1], linspace(0,2*tau,256)');
cmap=colormap(ax2, 'turbo');
caxis([0 thisDeltamax])
cb = colorbar(ax2);
cb.Ticks = [0, tau, 2*tau];
cb.TickLabels = {'0', '$\Delta=\tau$', '$\Delta=2\tau$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 18;
set(h, 'Visible','off')
axis xy
hold on

y_lasts = [];

for Dindx = inds

    Delta = Deltas(Dindx);
    this_phase_vars = phase_vars(:, Dindx);
    bool = this_phase_vars~=0;

    this_deltas = deltas(bool);
    this_y = this_phase_vars(bool);

    plot(ax2, [0 this_deltas], [this_y(1); this_y], Color = cmap(min([ceil(256*Delta/(thisDeltamax)), 256]), :), LineWidth = 1);
    hold(ax2, 'on');

    y_lasts = [y_lasts, this_y(end)];

end

% Connect end points
plot([0, Deltas(inds)], [0, y_lasts], '--', color = 'k')

% Marker at tau
[~, indx] = min(abs(Deltas(inds)-tau));
plot([tau, tau], [0, y_lasts(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')

% Phase variance saturation
yline(phase_var_sat, '--', color = 'k', LineWidth = 1, HandleVisibility='off')

% Axis limits
xlim([0 deltas(end-1)])
ylim([0 0.15])

% Axis ticks
yticks([0, phase_var_sat])
yticklabels({'', '$\frac{2}{5}\gamma^2 A^2 R^2$'})
xticks([0, tau])
xticklabels({'0', '$\tau$'})
set(ax2, 'TickLabelInterpreter', 'latex')
ax2.FontSize = 20;
ax2.YAxis.FontSize = 18;

xlabel('$\delta \rightarrow$', Position=[48,-0.002],Interpreter='latex', FontSize=20)
ylabel('$\langle \Delta\phi^2 \rangle \rightarrow$', Position=[-0.64 0.144] , Interpreter='latex', FontSize=18)


saveas(f2, 'sphere_phase_var_small_delta_1.png')

%% Plot phase variance + linear as function of delta 

% Predicted gradient
pred_grad = (-2/3)*((gamma*1e8)^2)*((A*1e-6)^2)*(2*1e-9)*(1e-3);

thisDeltamax = 56; % Maximum DELTA value to plot for
Imax = find(Deltas>thisDeltamax, 1, 'first');
inds = 1:Imax;

f3=figure;
f3.Position  = [286.6 97 916.8 657.6];
ax3 = axes;

% Create colorbar
h = imagesc(ax3, [0 2*tau], [0 1], linspace(0,2*tau,256)');
cmap=colormap(ax3, 'turbo');
caxis([0 thisDeltamax])
cb = colorbar(ax3);
cb.Ticks = [0, tau, 2*tau];
cb.TickLabels = {'0', '$\Delta=\tau$', '$\Delta=2\tau$'};
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 18;
set(h, 'Visible','off')
axis xy
hold on

y_lasts = [];

for Dindx = inds

    Delta = Deltas(Dindx);
    this_phase_vars = phase_vars(:, Dindx);
    bool = this_phase_vars~=0;

    this_deltas = deltas(bool);
    this_y = this_phase_vars(bool)-1*pred_grad*deltas(bool)';

    plot(ax3, [0 this_deltas], [this_y(1); this_y], Color = cmap(min([ceil(256*Delta/(thisDeltamax)), 256]), :), LineWidth = 1);
    hold(ax3, 'on');

    y_lasts = [y_lasts, this_y(end)];

end

% Connect end points
plot([0, Deltas(inds)], [0, y_lasts], '--', color = 'k')

% Marker at tau
[~, indx] = min(abs(Deltas(inds)-tau));
plot([tau, tau], [0, y_lasts(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')

% Phase variance saturation
yline(phase_var_sat, '--', color = 'k', LineWidth = 1, HandleVisibility='off')

% Axis limits
xlim([0 deltas(end-1)])
ylim([0 0.25])

% Axis ticks
yticks([0, phase_var_sat])
yticklabels({'', '$\frac{2}{5}\gamma^2 A^2 R^2$'})
xticks([0, tau])
xticklabels({'0', '$\tau$'})
set(ax3, 'TickLabelInterpreter', 'latex')
ax3.FontSize = 20;
ax3.YAxis.FontSize = 18;

xlabel('$\delta \rightarrow$', Position=[48,-0.003],Interpreter='latex', FontSize=20)
ylabel('$\langle \Delta\phi^2 \rangle + \frac{2}{3} \gamma^2 A^2 D \delta \rightarrow$', Position=[-0.64 0.204] , Interpreter='latex', FontSize=18)


saveas(f3, 'sphere_phase_var_small_delta_2.png')

%% Markers for sequences

R=10;

seq_deltas = ...[2,2,2,2,2,2,2,2,2,2]; 
            [3.5, 10.5, 23.5, 14, 20, 11, 18];

seq_Deltas = ...[15, 20, 30, 40, 50, 40, 60, 80, 100, 120]; 
            [26.9, 33.9, 46.9, 37.4, 43.4, 34.4, 41.4];

seq_bvals = ...[1000, 1250, 1500, 1750, 2000, 1000, 1250, 1500, 1750, 2000]; 
            [90, 500, 1500, 2000, 3000, 1000, 1800];


% seq_names = ...{
%     'Delta 15 b1000',...
%     'Delta 20 b1250',...
%     'Delta 30 b1500',...
%     'Delta 40 b1750',...
%     'Delta 50 b2000',...
%     'Delta 40 b1000',...
%     'Delta 60 b1250',...
%     'Delta 80 b1500',...
%     'Delta 100 b1750',...
%     'Delta 120 b2000',...
%     };

seq_names =  {'Classic b90',...
     'Classic b500',...
     'Classic b1500',...
     'Classic b2000',...
     'Classic b3000',...
     'Fast b1000',...
     'Fast b1800'};

cols = [colororder; colororder];

f4=figure;
f4.Position  = [286.6 97 916.8 657.6];
ax4 = axes;

for seqindx = 1:length(seq_deltas)

    delta = seq_deltas(seqindx);
    Delta = seq_Deltas(seqindx);

    thisG = A/delta;
    bval = stejskal(delta, Delta, G=thisG);

    % sphere phase variance
    sphere_signal = sphereGPD(delta, Delta, thisG, R, D);
    pv_sphere = -2*log(sphere_signal);

    % ball phase variance
    ball_signal = ball(bval, D);
    pv_ball = -2*log(ball_signal);   

    plot([Delta, Delta], [pv_ball, pv_sphere], '-', color = cols(seqindx,:), DisplayName = seq_names{seqindx})
    hold on
    scatter(Delta, pv_ball, '*', MarkerEdgeColor = cols(seqindx,:), DisplayName= 'Gaussian', HandleVisibility='off')
    scatter(Delta, pv_sphere, 'o', 'filled', ...
        MarkerFaceColor = cols(seqindx,:), MarkerEdgeColor = cols(seqindx,:), DisplayName= 'Sphere', HandleVisibility='off')
end


% Axis limits
xlim([0, thisDeltamax])
ylim([0 0.5])

% Connect start points of each curve
plot(deltas, phase_vars_firsts, '--', color = 'k', HandleVisibility='off');

% Mark tau and 2*tau
[~, indx] = min(abs(deltas-tau));
plot([tau, tau], [0, phase_vars_firsts(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')

[~, indx] = min(abs(deltas-2*tau));
plot([2*tau, 2*tau], [0, phase_vars_firsts(indx)], '--', LineWidth = 1, color = 'k', HandleVisibility='off')


% Phase variance saturation
yline(phase_var_sat, '--', color = 'k', LineWidth = 1, HandleVisibility='off')

% Phase variance linear growth (delta, Delta << tau)
plot(Deltas, alpha*(Deltas-deltas(1)/3), '--', color = 'k', LineWidth = 1, HandleVisibility='off')
text(26.5, 0.47, '$2 \gamma^2 A^2 D ( \Delta - \frac{\delta}{3})$', Interpreter='latex', FontSize=17)

% Axis ticks
yticks([0, phase_var_sat])
yticklabels({'', '$\frac{2}{5}\gamma^2 A^2 R^2$'})
xticks([0, tau, 2*tau])
xticklabels({'0', '$\tau$', '$2\tau$'})
set(ax4, 'TickLabelInterpreter', 'latex')
ax4.FontSize = 20;
ax4.YAxis.FontSize = 18;

legend(NumColumns=1, FontSize=14, Location="northwest", Interpreter='latex')

grid on

xlabel('$\Delta \rightarrow$', Position=[56,-0.0065],Interpreter='latex', FontSize=20)
ylabel('$\langle \Delta\phi^2 \rangle \rightarrow$', Position=[-0.64 0.484] , Interpreter='latex', FontSize=18)

saveas(f4, 'phase_var_sequences.png')



% === Look at signals

f5 = figure;
f5.Position  = [286.6 97 916.8 657.6];
ax5 = axes;


for seqindx = 1:length(seq_deltas)

    delta = seq_deltas(seqindx);
    Delta = seq_Deltas(seqindx);
    bval = seq_bvals(seqindx);

    thisG = stejskal(delta, Delta, bval=bval);

    % sphere phase variance
    sphere_signal = sphereGPD(delta, Delta, thisG, R, D);

    % Apparent diffusivity
    sphere_Dapp = (-1/bval)*log(sphere_signal)*1e3;

    % ball phase variance
    ball_signal = ball(bval, D);
    ball_Dapp = D;

    plot([Delta, Delta], [ball_signal, sphere_signal], '-', color = cols(seqindx,:), DisplayName = seq_names{seqindx})
    hold on
    scatter(Delta, ball_signal, '*', MarkerEdgeColor = cols(seqindx,:), DisplayName= 'Gaussian', HandleVisibility='off')
    scatter(Delta, sphere_signal, 'o', 'filled', ...
        MarkerFaceColor = cols(seqindx,:), MarkerEdgeColor = cols(seqindx,:), DisplayName= 'Sphere', HandleVisibility='off')

end

xlim([3*tau/4, (2.25)*tau])
xticks([0, tau, 2*tau])
xticklabels({'0', '$\tau$', '$2\tau$'})

ylim([-0.05 1])
yticks([0 0.25 0.5 0.75 1])
ylabel('Signal', Interpreter='latex')
set(ax5, 'TickLabelInterpreter', 'latex')
ax5.FontSize = 20;
ax5.YAxis.FontSize = 18;

grid on

legend(Location='southwest', FontSize=14, Interpreter='latex')

xlabel('$\Delta \rightarrow$', Position=[56,-0.065],Interpreter='latex', FontSize=20)

saveas(f5, 'signal_sequences.png')




% == Look at apparent diffusivity of sphere compartment


f6 = figure;
f6.Position  = [286.6 97 916.8 657.6];
ax6 = axes;


for seqindx = 1:length(seq_deltas)

    delta = seq_deltas(seqindx);
    Delta = seq_Deltas(seqindx);
    bval = seq_bvals(seqindx);

    thisG = stejskal(delta, Delta, bval=bval);

    % sphere phase variance
    sphere_signal = sphereGPD(delta, Delta, thisG, R, D);

    % Apparent diffusivity
    sphere_Dapp = (-1/bval)*log(sphere_signal)*1e3;

    scatter(Delta, sphere_Dapp, 'o', 'filled', ...
        MarkerFaceColor = cols(seqindx,:), MarkerEdgeColor = cols(seqindx,:), DisplayName = seq_names{seqindx})
    hold on
    plot([0 Delta], [sphere_Dapp sphere_Dapp], '--', color = cols(seqindx,:), HandleVisibility='off')
end

xlim([3*tau/4, (2.25)*tau])
xticks([0, tau, 2*tau])
xticklabels({'0', '$\tau$', '$2\tau$'})

ylim([-0.0 2])
yticks([0 0.2 0.4 0.6 0.8 1 1.5 2])

ylabel('Sphere apparent diffusivity ($\times 10^{-3}$ mm$^2$/s)', Interpreter='latex')
set(ax6, 'TickLabelInterpreter', 'latex')
ax6.FontSize = 20;
ax6.YAxis.FontSize = 16;

grid on

legend(Location='northeast', FontSize=14, Interpreter='latex')

xlabel('$\Delta \rightarrow$', Position=[56,-0.03],Interpreter='latex', FontSize=20)

saveas(f6, 'Dapp_sequences.png')