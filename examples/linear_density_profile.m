%% GRAVITY COEFFICIENTS OF LINEAR DENSITY PLANET
% Example and test of the CMSPlanet class. We construct and converge a
% model of a rotating fluid planet with a density profile that is linear in the
% mean radius of a level surface.
%
% The comparison is with a 5th-order Zharkov and Trubitsyn model. I take the
% numbers from Table 1 of Hubbard (2013). A complication is that the Z+T theory
% works in powers of the small parameter m rather than q. A CMS-based model
% would have to iteratively converge on a the desired values of m and s_i. But I
% will use the value of q already found by Hubbard.

%% Prepare workspace
clear
clc
close all
si = setFUnits;
G = si.gravity;

%% Set up a CMS object and give it a density profile linear in lambda to start
N = 128;
nx = 128;
cmp = CMSPlanet;
cmp.opts.xlayers = nx;

% The lambda array follows Hubbard's IDL code
dl = 1/(N - 1);
lam(1) = 1;
lam(2) = 1 - dl/2;
for k=3:N
    lam(k) = lam(k-1) - dl;
end
cmp.ai = lam;

% The deltas array follows Hubard's IDL code
deltas(N) = 0;
deltas(2:end) = dl;
cmp.rhoi = cumsum(deltas);

cmp.qrot = 0.088822426; % Hubbard (2013) Table 1

%% Relax to hydrostatic equilibrium
cmp.relax_to_HE;

%% After initial relaxation, we iteratively fix deltas ss and re-relax
for iter=1:20
    fprintf('\n  Fixing density profile to mean radii - iteration %i\n', iter)

    % Following setup_linear_iterated.pro
    sos0 = cmp.si/cmp.s0;
    dsos0 = [-diff(sos0); sos0(end)];
    dlambda = [-diff(cmp.CMS.lambdas); cmp.CMS.lambdas(end)];
    delta0 = [0; ones(N-1,1)/(N-1)];
    new_deltas = delta0.*dsos0./dlambda;
    
    delta_deltas = sum(abs(new_deltas - cmp.CMS.deltas));
    fprintf('  deltas array modified by < %g.\n\n',delta_deltas)
    cmp.rhoi = cumsum(new_deltas);
    if delta_deltas < 1e-12, break, end
    cmp.relax_to_HE;
end
fprintf('  Density profile fixed.\n')

%% Compare computed and analytic density structure
q = cmp.qrot;
m = cmp.mrot;

% Zharkov & Trubistyn (1978) Table 3.1
ZT5 = [nan; 0.0830; 1.4798; 5.929; 3.497; 2.52; 2.4; nan; nan];
% Hubbard (2013) Table 1
H13_128 = [q; 0.082999915; 1.4798138; 5.9269129; 3.4935680; 2.5493209;...
    2.1308951; 1.9564143; 1.9237724];

% CMSPlanet
CMP = [q; m; cmp.J2*1e2; -cmp.J4*1e4; cmp.J6*1e5; -cmp.J8*1e6;...
    cmp.Js(6)*1e7; -cmp.Js(7)*1e8; cmp.Js(8)*1e9];

% Make it a table
T = table(ZT5, H13_128, CMP);
T.Properties.RowNames = {'q','m','J2x10^2','-J4x10^4','J6x10^5','-J8x10^6',...
    'J10x10^7','-J12x10^8','J14x10^9'};

% Display
format long
format compact
fprintf('\n')
disp(T)
format

%% Save and deliver
save('linear_density_model.mat', 'cmp', 'T')
