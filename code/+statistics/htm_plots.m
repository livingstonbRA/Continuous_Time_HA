
kappa1 = kappa1(:);
kappa2 = kappa2(:);
phtm = phtm(:);
whtm = whtm(:);
ra = ra(:);
mpc = mpc(:);

kappa2s = [0.1, 0.5, 1, 2];

variables = {phtm, whtm, ra, mpc};
labels = {'PHtM', 'WHtM', 'r_a', 'QMPC'};

for iv = 1:4
    var = variables{iv};
    label = labels{iv};

    for k2 = kappa2s
        scatter(kappa1(kappa2==k2), var(kappa2==k2))
        hold on
    end
    hold off

    xlabel('kappa1')
    ylabel(label)
    legend('kappa2 = 0.1', 'kappa2 = 0.5', 'kappa2 = 1', 'kappa2 = 2')
    set(gcf,'color','w');
    fname = sprintf('output/%s_plot.png', label);
    saveas(gcf, fname)

end