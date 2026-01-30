% Make example graphs for restricted diffusion


tau = 40;
ts = 0.0:0.2:200;

FS = 22;

f=figure;
f.Position = [384.2 158.6 766.4 545.6];

% plot(ts, 1./(1+2*(ts/tau).^1), LineWidth = 1.2)
plot(ts, func(ts, tau), LineWidth = 1.2)
hold on
scatter(0, 1, 10,  'filled', MarkerFaceColor = [0, 0.4470, 0.7410])
plot([1.5*tau 1.5*tau], [0, func(1.5*tau, tau)], '--', color = 'k')


text(14, 0.95, '$\approx \frac{R^2}{5} - \frac{Dt}{3}$', Interpreter='latex', FontSize=FS)
text(155, 0.25, '$\sim \frac{R^4}{Dt}$', Interpreter='latex', FontSize=FS)


ylim([0 1.5])
xlim([0 max(ts)])

xticks([0 1.5*tau])
xticklabels(["{0}", "${\tau}$"])
yticks([1])
yticklabels(["${\frac{R^2}{5}}$"])

ax = gca();
ax.FontSize = FS;
set(ax, 'TickLabelInterpreter', 'latex')

xlabel('$t \rightarrow$', Position=[max(ts)-6,-0.01],Interpreter='latex', FontSize=18)
ylabel('$\langle \bar{x}^2 \rangle \rightarrow$', Position=[-2 1.38] , Interpreter='latex', FontSize=FS)

box on


saveas(f, 'time-averaged position variance.png')



function vals = func(ts, tau)

vals = exp(-0.5*ts/tau);


end