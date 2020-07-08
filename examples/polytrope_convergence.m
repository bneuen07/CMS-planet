function T = polytrope_convergence(N,nx,zgrid,testname)

if nargin < 4
    fprintf('Usage:\n\tpolytrope_convergence(N,nx,zgrid,testname)\n')
    fprintf('\tzgrid(n) returns a lambda grid.\n')
    return
end

% Create a table to store experiment results
vars = {'N','nx','J2','runtime'};
tps = {'double','double','double','double'};
T = table('Size', [length(N), length(vars)],...
    'VariableTypes', tps,...
    'VariableNames', vars);

% Run cms, this may take a while
fprintf("Running CMS, this may take a while!\n")

for j=1:length(N)
    fprintf('working on row %d of %d (N=%d)...',j,height(T),N(j));
    T.N(j) = N(j);
    T.nx(j) = nx(j);
    zvec = zgrid(N(j));
    
    % Set up a CMS Planet (mass and radius arbitrary)
    M = 317.8*5.97e24;
    R = 71492*1e3;

    cmp = CMSPlanet;
    cmp.ai = R*zvec;
    cmp.rhoi = -zvec.^2 + 1;
    cmp.mass = M;
    cmp.radius = R;
    cmp.renormalize();
    cmp.qrot = 0.089195487; % Hubbard (2013) Table 5
    cmp.P0 = 0.0;

    % Construct a polytrope of index 1 to represent the planet's eos
    n = 1;
    K = 3.4e5; % K has no effect on the Js
    eos = barotropes.Polytrope(K, n);
    cmp.eos = eos;

    %% Relax to desired barotrope
    cmp.opts.MaxIterBar = 40;
    cmp.opts.xlayers = nx(j);
    cmp.opts.dJtol = 1e-8;
    cmp.opts.verbosity = 0;
    tic
    cmp.relax_to_barotrope;
    T.runtime(j) = toc;
    T.J2(j) = cmp.J2;
    fprintf('done (%d sec).\n',T.runtime(j))
    save(testname, 'T')
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));
save(testname, 'T')

