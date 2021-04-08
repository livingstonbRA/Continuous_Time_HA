function scf = scf2019struct()
    scf = struct();
    
    scf.quarterly_earnings = 79181 / 4;
    scf.median_totw = 1.54;
    scf.median_liqw = 0.05;
    scf.htm = 0.39;
    scf.phtm = 0.135;
    scf.mean_totw = 3.5;

    %% Notes
    % HtM defined as household with b < y / 6
    % PHtM defined as household with (a + b) < y / 6
    
end