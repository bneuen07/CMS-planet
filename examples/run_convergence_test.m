%% How many density layers and shape layers are required for desired precision?
clear
clc
close all

%% Choose layer numbers to investigate
N = [128];
nx = [32];

%% The experiment MAY be mildly sensitive to the chosen density model
rho0 = 0.199;
p=[-9.1352277e+02 -2.0616194e+03 -2.5768977e+02 -6.8877550e+01 8.6817818e+03 4.2076235e+03 -1.3579913e+04 0.0 3.9924165e+03];
p=[  -2.0620547e+03 1.9490798e+03 1.7669895e+03 -6.5775595e+03 1.1152762e+04 -4.0190954e+03 -6.5570643e+03 0.0 4.3471119e+03];
cmp = generators.polynomial(4096, p);
z4k = cmp.ai/cmp.a0;
d4k = cmp.rhoi;
zstrat = @(n)lambdas.MWH19(n, z4k, d4k);
zgridname = 'MWH19';
mdl = @(n)generators.polynomial(n,p,zstrat);
qrot = 0.1574;

%% Create a table to store experiment results
vars = {'N','nx','J2','runtime','cmp'};
tps = {'double','double','double','double','CMSPlanet'};
T = table('Size', [length(N)*length(nx), length(vars)],...
    'VariableTypes', tps,...
    'VariableNames', vars);
[~,fname,~] = fileparts(tempname);

%% Relax CMPs, this may take a while
fprintf("Relaxing CMPs, this may take a while!\n")

for j=1:length(N)
    for k=1:length(nx)
        row = (j-1)*length(nx) + k;
        fprintf('working on row %d of %d (N=%d, nx=%d)...',row,height(T),...
            N(j),nx(k));
        T.N(row) = N(j);
        T.nx(row) = nx(k);
        cmp = mdl(N(j));
        cmp.opts.xlayers = nx(k);
        cmp.qrot = qrot;
        trun = cmp.relax_to_HE;
        T.runtime(row) = trun;
        T.J2(row) = cmp.J2;
        T.cmp(row) = cmp;
        fprintf('done.\n')
    end
end
fprintf('All done. (%s)\n',seconds2human(sum(T.runtime)));
save(zgridname, 'T')
%% Examine the convergence behavior
