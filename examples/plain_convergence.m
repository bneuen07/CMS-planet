%% How many layers before convergence?
clear
clc
close all

%% Choose layer numbers to investigate
N = [512, 1024];
testname = 'MWH19_nint1';

%% The experiment MAY be mildly sensitive to the chosen density model
p=[  -2.0620547e+03 1.9490798e+03 1.7669895e+03 -6.5775595e+03 1.1152762e+04 -4.0190954e+03 -6.5570643e+03 0.0 4.3471119e+03];
cmp = generators.polynomial(4096, p);
z4k = cmp.ai/cmp.a0;
d4k = cmp.rhoi;
zstrat = @(n)lambdas.MWH19(n, z4k, d4k);
mdl = @(n)generators.polynomial(n,p,zstrat);
qrot = 0.08;

%% Create a table to store experiment results
vars = {'N','J2','runtime'};
tps = {'double','double','double'};
T = table('Size', [length(N), length(vars)],...
    'VariableTypes', tps,...
    'VariableNames', vars);

%% Run cms, this may take a while
fprintf("Running CMS, this may take a while!\n")

for j=1:length(N)
    row = j;
    fprintf('working on row %d of %d (N=%d)...',row,height(T),N(j));
    T.N(row) = N(j);
    cmp = mdl(N(j));
    cmp.opts.xlayers = N(j);
    cmp.qrot = qrot;
    trun = cmp.relax_to_HE;
    T.runtime(row) = trun;
    T.J2(row) = cmp.J2;
    fprintf('done.\n')
    save
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));
save(testname, 'T')
