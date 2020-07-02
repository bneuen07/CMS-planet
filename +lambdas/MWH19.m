function lams = MWH19(N, zvec, dvec)
%MWH19 Return a lambda distribution following prescription of Militzer et al 2019.
%    lambdas = MWH19(N, zvec, dvec) returns an N-vector of normalized radii
%    following the optimized grid idea from Militzer et al. (2019). This lambda
%    distribution is model-depdendent, which is why we need an idea of the density
%    profile to begin with, supplied in zvec and dvec. The prescription calls for
%    constructing a geometric sequence of density values rho(i) between dvec(1)
%    and dvec(end) such that rho(i+1)/rho(i) is constant. Then the supposedly
%    optimal lambda grid is the one where dvec at lambda_i equals rho(i).

if nargin < 3
    fprintf('Usage:\n\tMWH19(N, zvec, dvec)\n')
    fprintf('\tMake sure zvec and dvec are monotonic indexed from surface in.\n')
    return
end

f = (dvec(end)/dvec(1))^(1/(N - 1));
assert(f > 1) % just making sure
ros = ones(N,1);
ros(1) = dvec(1);
for k=2:N
    ros(k) = ros(k-1)*f;
end
ros(N) = dvec(end); % it's already "equal" but interpolation needs identical.

lams = interp1(dvec, zvec, ros);
