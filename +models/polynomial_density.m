function cmp = polynomial_density(N, x, lamstrat, forcezero)
%POLYNOMIAL_DENSITY Planet model with a polynomial density structure.
%    POLYNOMIAL_DENSITY(N, x) returns an N-layer CMSPlanet object with a density
%    structure given by a polynomial whose coefficients are in x. The density will
%    be:
%
%    rho(lamb) = x(1) + x(2)*lamb + x(3)*lamb^2 + x(4)*lamb^3 + ...
%
%    The degree of the polynomial is selected to agree with the length of the
%    vector x given in input. Minimum length is one (constant density). The
%    central density, x(1), must be positive. To have zero density at the top
%    layer, the sum of the coeffients must be zero (a warning will appear if it is
%    not so). Note that the actual rhoi quantity of output cmp model will be
%    scaled to match the mass and radius of the model.
%
%    POLYNOMIAL_DENSITY(N, x, lamstrat) where lamstrat is a 2-element vector lets
%    you specify the layer spacing strategy. Approximately lamstrat(1) of
%    available layers will be distributed in the top lamstrat(2) of the planet.
%    For example, passing [3/4, 0.2] concentrates the layers heavily in the top
%    20% of the planet, leaving about N/4 layers to fill the bottom 80%. A single
%    half-width layer of zero density is still reserved for the surface. To use
%    the default spacing pass lambda=[].
%
%    POLYNOMIAL_DENSITY(N, x, lamstrat) where lamstrat is a function handle lets
%    you specify the lambda spacing completely. Pass a handle to a function that
%    takes a single scalar integer (number of layers) and returns a vector of that
%    length with values in the interval (0, 1], for normalized layer radii. For
%    example, to set layers with equally spaced radii use
%    lamstrat=@(n)linspace(1,1/n,n).
%
%    POLYNOMIAL_DENSITY(..., forcezero) if forcezero=true forces the top layer to
%    have zero density.

narginchk(2,4)
if ((nargin == 2) || isempty(lamstrat)), lamstrat = [2/3, 1/2]; end
if ((nargin < 4) || isempty(forcezero)), forcezero = false; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'vector',}, '', 'x', 2)
assert(x(1) > 0, 'The central density x(1) must be positive.')
validateattributes(lamstrat, {'numeric','function_handle'}, {}, '', 'lamstrat', 3)
if isnumeric(lamstrat)
    validateattributes(lamstrat, {'numeric'},...
        {'vector', 'numel', 2, '>', 0, '<', 1}, '', 'lamstrat', 3)
end
validateattributes(forcezero, {'logical'}, {'scalar'}, '', 'forcematch', 4)

cmp = CMSPlanet(N);

if (isa(lamstrat, 'function_handle'))
    lambdas = lamstrat(N);
    assert(isnumeric(lambdas) && isvector(lambdas) && (numel(lambdas) == N),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    assert(all(lambdas > 0) && all(lambdas <= 1),...
        '@lamstrat(N) must return a vector of length N with values in (0,1].')
    cmp.cms.lambdas = lambdas;
else
    n1 = fix(lamstrat(1)*(N - 1));
    n2 = N - n1 - 1;
    dl1 = lamstrat(2)/(n1 - 1);
    dl2 = (1 - lamstrat(2))/(n2 + 1);
    lam1 = linspace(1 - dl1/2, (1 - lamstrat(2)), n1);
    lam2 = linspace((1 - lamstrat(2)) - dl2, dl2, n2);
    cmp.cms.lambdas = [1, lam1, lam2]';
end

degree = length(x) - 1;

rhoi = repmat(x(1), N, 1);
for n = 1:degree;
    rhoi = rhoi + x(n+1)*power(cmp.cms.lambdas, n);
end

if forcezero; rhoi(1) = 0.0; end

cmp.cms.deltas = [rhoi(1); diff(rhoi)];

% Some warning signs
if any(rhoi < 0)
    warning('Some densities are negative!')
end
if abs(sum(x)) > eps
    warning('Density in the top layer is not be zero!')
end
if any(cmp.cms.deltas(2:end) < 0)
    warning('The density structure is not monotonic!')
end

end
