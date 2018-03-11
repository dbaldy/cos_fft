% Compute price for European vanilla calls and puts with the Cos-FFT method described in Fang & Oosterlee (2008).
%
% danilo.zocco@gmail.com, 2018-03-11

% Option parameters
S0 = 100;   % Stock price
r = 0;      % Risk-free interest rate
mu = r;     % Drift (should be r-q, with q being a constant yield)
T = 0.1;    % Maturity in years
K = 80:5:120;          % Vector of strike prices
option_type = 'call';  % Call or Put (case insensitive)

% Heston parameters
u0 = 0.0175;      % Current variance
u_bar = 0.0398;   % Long-run variance
lambda = 1.5768;  % Rate of mean reversion
eta = 0.5751;     % Vol of vol
rho = -0.5711;    % Correlation of Wiener Processes

% Cos-FFT parameters
N = 160;
k = 0:N-1;
w = [0.5, ones(1,N-1)];
x = log(S0./K);

% Feller condition (Does not affect the model but gets mentioned in Fang & Oosterlee (2008))
% if 2*u_bar*lambda > eta^2
%     disp('Feller condition is respected')
% else
%     disp('Feller condition is not respected')
% end

% Truncation range
    % Cumulants c1 and c2
    c1 = mu*T + (1-exp(-lambda*T))*(u_bar-u0)/(2*lambda) - 0.5*u_bar*T;

    c2 = (eta*T*lambda*exp(-lambda*T)*(u0-u_bar)*(8*lambda*rho-4*eta)...
        + lambda*rho*eta*(1-exp(-lambda*T))*(16*u_bar-8*u0)...
        + 2*u_bar*lambda*T*(-4*lambda*rho*eta+eta^2+4*lambda^2)...
        + eta^2*((u_bar-2*u0)*exp(-2*lambda*T)+u_bar*(6*exp(-lambda*T)-7)+2*u0)...
        + 8*lambda^2*(u0-u_bar)*(1-exp(-lambda*T)))/(8*lambda^3);

    % Limits of support a and b
        % Might need to add x to a and b, and use L=13
    a = c1-12*sqrt(abs(c2));
    b = c1+12*sqrt(abs(c2));
    
% Characteristic function
D = @(omega) sqrt((lambda-1i*rho*eta*omega).^2 + (omega.^2+1i*omega)*eta^2);
G = @(omega) (lambda-1i*rho*eta*omega-D(omega))./(lambda-1i*rho*eta*omega+D(omega));

CharacteristicFunction = @(omega) ...
    exp(1i*omega*mu*T+u0*((1-exp(-D(omega)*T))./(1-G(omega).*exp(-D(omega)*T))).*(lambda-1i*rho*eta*omega-D(omega))/eta^2) ...
    .* exp((T*(lambda-1i*rho*eta*omega-D(omega))-2*log((1-G(omega).*exp(-D(omega)*T))./(1-G(omega))))*lambda*u_bar/eta^2);

% Parameter Omega
omega = k*pi/(b-a);

% U_k
Chi = @(argin) (cos(omega*(argin(2)-a))*exp(argin(2))-cos(omega*(argin(1)-a))*exp(argin(1)) ...
        + omega.*sin(omega*(argin(2)-a))*exp(argin(2))-omega.*sin(omega*(argin(1)-a))*exp(argin(1))) ...
        ./(1+omega.^2);

Psi = @(argin) [argin(2)-argin(1), ...      % For k = 0
        (sin(omega(2:end)*(argin(2)-a))-sin(omega(2:end)*(argin(1)-a)))./omega(2:end)];   % For k > 0

if strcmpi(option_type,'call')
    argin = [0,b];
    U_k = 2*(Chi(argin)-Psi(argin))/(b-a);
elseif strcmpi(option_type,'put')
    argin = [a,0];
    U_k = 2*(-Chi(argin)+Psi(argin))/(b-a);
else
    disp('This option type is not supported.')
end


% Options prices
v = K * exp(-r*T) .* real(w.*CharacteristicFunction(omega).*U_k*exp((x-a).'*1i*omega).');

disp([K;v]) % Display a matrix with strikes in top row and option prices in bottom row.







