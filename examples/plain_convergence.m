%% How many layers before convergence?
clear
clc
close all

%% Choose layer numbers and lambda grid to investigate
N = [16,32];
nx = N;
zgrid = @(n)linspace(1, 1/n, n); % equally spaced radii
testname = 'linear_allx';

%% This experiment may be moderately sensitive to the density profile used
dprof = @(z)-z.^2 + 1; % simple planetary quadratic

%% Create a table to store experiment results
vars = {'N','nx','J2','runtime'};
tps = {'double','double','double','double'};
T = table('Size', [length(N), length(vars)],...
    'VariableTypes', tps,...
    'VariableNames', vars);

%% Run cms, this may take a while
fprintf("Running CMS, this may take a while!\n")

for j=1:length(N)
    fprintf('working on row %d of %d (N=%d)...',j,height(T),N(j));
    T.N(j) = N(j);
    T.nx(j) = nx(j);
    
    zvec = zgrid(N(j));
    dvec = dprof(zvec);
    qrot = 0.1;
    tic; [Js, out] = cms(zvec, dvec, qrot, 'xlayers', nx(j), 'tol', 1e-10);
    T.runtime(j) = toc;
    T.J2(j) = Js(2);
    fprintf('done (%d sec).\n',T.runtime(j))
    save
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));
save(testname, 'T')
