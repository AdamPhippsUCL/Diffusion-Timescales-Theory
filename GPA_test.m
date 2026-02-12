% Script to test validity of GPA for delta<<tau<Delta



R=10*10^(-6);
gamma=2.675*10^8; %



f=figure;
f.Position  = [680   150   759   588];

%%


ss = -R:(R/200):R;


As = (1:1:240)*10^(-3)*(20*10^-3);  % G*dela=ta

sin_exps = zeros(size(As));
cos_exps = zeros(size(As));
var_phis = zeros(size(As));

sphere_signals = zeros(size(As));

for Aindx = 1:length(As) % 100 mT/m * 10 ms

    A = As(Aindx);

    % % Test distribution shape
    % figure
    % for s = -2*R:0.2:2*R
    %     scatter(s, pdf(s,R))
    %     hold on
    % 
    % end



    % Distribution sum
    dist_sum = 0;
    for sindx = 1:length(ss)
        s = ss(sindx);
        dist_sum = dist_sum + pdf(s, R);
    
    end


    
    
    % Distribution variance <phi^2>
    
    var_s = 0;
    
    for sindx = 1:length(ss)
        
        s = ss(sindx);
    
        var_s = var_s + (s^2)*pdf(s, R)/dist_sum;
    
    end
    
    std_s = sqrt(var_s);
    
    var_phi = (gamma^2*A^2)*var_s;
    
    var_phis(Aindx) = var_phi;
    
    
    % Cosine expectation
    
    cos_exp = 0;
    sin_exp = 0;

    for sindx = 1:length(ss)
        
        s = ss(sindx);
    
        cos_exp = cos_exp + cos(gamma*A*s)*pdf(s, R)/dist_sum;
        sin_exp = sin_exp + sin(gamma*A*s)*pdf(s, R)/dist_sum;
    
    end
    
    cos_exps(Aindx) = cos_exp;
    sin_exps(Aindx) = sin_exp;   



   % sphere_signals(Aindx) = sphereGPD(0.1, 200, 100*A/0.1, R, 2);

end


exp_exps = sqrt(cos_exps.^2+sin_exps.^2);

% plot(As, var_phis)
plot(gamma*As*R, exp(-0.5*var_phis), '--', DisplayName='$e^{\frac{-\langle \Delta \phi^2 \rangle}{2}}$', LineWidth = 1.8)
hold on
plot(gamma*As*R, abs(exp_exps), DisplayName = '$\langle e^{-i\Delta\phi} \rangle$', LineWidth = 1.8)
% plot(As, sphere_signals)
legend(Interpreter='latex', FontSize=20)

xlabel('$\gamma A R$ (radians)', Interpreter= 'latex')
ylabel('{Attenutation factor}', Interpreter= 'latex')

ax = gca();
ax.FontSize = 14;

xlim([0, 8])
ylim([0, 1.1])


saveas(f, 'GPA_short_pulse_long_separation.png')




% % Error
% f2 = figure;
% plot(gamma*As*R, exp(-0.5*var_phis) - abs(exp_exps))
% 
% 
% xlabel('$\gamma A R$ (radians)', Interpreter= 'latex')
% ylabel('{Error}', Interpreter= 'latex')
% 
% ax = gca();
% ax.FontSize = 14;
% 
% xlim([0, 8])
% ylim([0, 1.1])


function pd = pdf(s, R)
    pd = (3/(160*R^6))*( 32*R^5 - 40*(R^3)*s^2 + 20*(R^2)*abs(s)^3 - abs(s)^5);
end


