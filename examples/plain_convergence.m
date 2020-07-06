function T = plain_convergence(N,nx,dprof,zgrid,testname)

if nargin < 5
    fprintf('Usage:\n\tplain_convergence(N,nx,dprof,zgrid,testname)\n')
    fprintf('\tzgrid(n) returns a lambda grid.\n')
    fprintf('\tdprof(z) returns a density profile.\n')
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
    dvec = dprof(zvec);
    qrot = 0.1;
    tic; Js = cms(zvec, dvec, qrot, 'xlayers', nx(j), 'tol', 1e-10);
    T.runtime(j) = toc;
    T.J2(j) = Js(2);
    fprintf('done (%d sec).\n',T.runtime(j))
    save(testname, 'T')
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));
save(testname, 'T')
