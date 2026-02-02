% Make figure for Gaussian diffusion


D=2;
times = [0, 1, 2, 4, 6, 8, 10];

xmax = 20;
xs = -xmax:xmax/100:xmax;
dx = xs(2)-xs(1);

f=figure;
ax = axes;

for tindx = 1:length(times)

    t = times(tindx);

    if t==0
        xline(0, DisplayName=['t=' num2str(t) 'ms'], LineWidth=1.1);
        hold on
        continue
    end

    sigma = sqrt(2*D*t);

    dist = normpdf(xs, 0, sigma);

    plot(xs, dx*dist, DisplayName=['t=' num2str(t) 'ms'], LineWidth=1.1);

end

legend(Interpreter='latex', FontSize=11);


ylim([0 0.05])
xlim([-xmax, xmax])
title(['$D =$ ' num2str(D) ' $\mu$m$^2$/ms'], Interpreter='latex', FontSize=14)
xlabel('$x$ ($\mu$m)', Interpreter='latex', FontSize=14)
yticks([0:0.01:0.05])
ylabel(['$P(x)$ ($\mu$m$^{-1}$)' ], Interpreter='latex', FontSize=14)

saveas(f, 'Gaussian_diffusion_profiles.png')


