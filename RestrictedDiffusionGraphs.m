% Make example graphs for restricted diffusion


tau = 40;
ts = 0.0:0.2:250;

FS = 23;

cols = colororder;

f=figure;
f.Position = [384.2 158.6 800 545.6];


plot(ts(1:end), func(ts(1:end), tau), LineWidth = 1.2, color = cols(1,:))
hold on
% plot(ts(201:450), func(ts(201:450), tau), '--', LineWidth = 1.2, color = cols(1,:))
% hold on
% plot(ts(451:end), func(ts(451:end), tau),LineWidth = 1.2, color = cols(1,:))

scatter(0, 1, 10,  'filled', MarkerFaceColor = [0, 0.4470, 0.7410])
plot([1.5*tau 1.5*tau], [0, func(1.5*tau, tau)], '--', color = 'k')


text(14, 0.95, '$\approx \frac{R^2}{5} - \frac{D  \Delta t}{3}$', Interpreter='latex', FontSize=FS)
text(180, 0.2, '$\sim \frac{R^4}{D} \frac{1}{\Delta t}$', Interpreter='latex', FontSize=FS)


ylim([0 1.5])
xlim([0 max(ts)])

xticks([0 1.5*tau])
xticklabels(["{0}", "${\tau}$"])
yticks([1])
yticklabels(["${\frac{R^2}{5}}$"])

ax = gca();
ax.FontSize = FS+1;
set(ax, 'TickLabelInterpreter', 'latex')

xlabel('$\Delta t \rightarrow$', Position=[max(ts)-6,-0.02],Interpreter='latex', FontSize=FS-2)
ylabel('$\langle \bar{x}^2 \rangle \rightarrow$', Position=[-2 1.38] , Interpreter='latex', FontSize=FS)

box on


saveas(f, 'time-averaged position variance.png')



function vals = func(ts, tau)

vals = exp(-0.5*ts/tau);


end