function y = TR2(Lift,span,Uinf,beta)

rho = 1.154;
le = span / 2;	% length from root to tip
N = 160;	    % Partition number
delta_S = linspace(le / N / 2, le / N / 2, N);	% half of partition, Equal distance

y = linspace(delta_S(1), le - delta_S(N), N);
z = linspace(0, 0, N);
fai = linspace(0, 0, N);

% ---- Geometric variable number ----
% Variables required to seek the Q.
for i = 1:N
	for j = 1:N
		ydash(i,j)  =   (y(i) - y(j)) .* cos(fai(j)) + (z(i) - z(j)) .* sin(fai(j));
		zdash(i,j)  = - (y(i) - y(j)) .* sin(fai(j)) + (z(i) - z(j)) .* cos(fai(j));
		y2dash(i,j) =   (y(i) + y(j)) .* cos(fai(j)) - (z(i) - z(j)) .* sin(fai(j));
		z2dash(i,j) =   (y(i) + y(j)) .* sin(fai(j)) + (z(i) - z(j)) .* cos(fai(j));

		R2puls(i,j)  = (ydash(i,j) - delta_S(j)).^2 + zdash(i,j).^2;
		R2minus(i,j) = (ydash(i,j) + delta_S(j)).^2 + zdash(i,j).^2;
		Rdash2puls(i,j)  = (y2dash(i,j) + delta_S(j)).^2 + z2dash(i,j).^2;
		Rdash2minus(i,j) = (y2dash(i,j) - delta_S(j)).^2 + z2dash(i,j).^2;
	end
end

for i = 1:N
	for j = 1:N
		Q(i,j) = - 1 / (2 * pi) *(((ydash(i,j) - delta_S(j)) / R2puls(i,j)- (ydash(i,j) + delta_S(j)) / R2minus(i,j)) * cos(fai(i) - fai(j))+ (zdash(i,j) / R2puls(i,j) - zdash(i,j) / R2minus(i,j)) * sin(fai(i) - fai(j))+ ((y2dash(i,j) - delta_S(j)) / Rdash2minus(i,j)- (y2dash(i,j) + delta_S(j)) / Rdash2puls(i,j)) * cos(fai(i) + fai(j))+ (z2dash(i,j) / Rdash2minus(i,j) - z2dash(i,j) / Rdash2puls(i,j))* sin(fai(i) + fai(j)));
	end
end

% ---- Normalization ----
% Variables required to seek q.
delta_sigma = delta_S ./ le;
eta = y ./ le;
etadash = ydash ./ le;
eta2dash = y2dash ./ le;
zeta = z ./ le;
zetadash = zdash ./ le;
zeta2dash = z2dash ./ le;
gamma2puls = R2puls / (le^2);
gamma2minus = R2minus / (le^2);
gammadash2puls = Rdash2puls / (le^2);
gammadash2minus = Rdash2minus / (le^2);
for i = 1:N
	for j = 1:N
		q(i,j) = - 1 / (2 * pi) *(((etadash(i,j) - delta_sigma(j)) / gamma2puls(i,j)- (etadash(i,j) + delta_sigma(j)) / gamma2minus(i,j)) * cos(fai(i) - fai(j))+ (zetadash(i,j) / gamma2puls(i,j) - zetadash(i,j) / gamma2minus(i,j))* sin(fai(i) - fai(j))+ ((eta2dash(i,j) - delta_sigma(j)) / gammadash2minus(i,j)- (eta2dash(i,j) + delta_sigma(j)) / gammadash2puls(i,j)) * cos(fai(i) + fai(j))+ (zeta2dash(i,j) / gammadash2minus(i,j) - zeta2dash(i,j) / gammadash2puls(i,j))* sin(fai(i) + fai(j)));
	end
end

% ---- elliptic loading aerodynamic force ----
% Vn is Induced vertical velocisy.
% Vn is constant when elliptical circulation distribution.
BendingMoment_elpl = 2 / 3 / pi * le * Lift;
InducedDrag_elpl = Lift^2 / (2 * pi * rho * Uinf^2 * le^2);
Vn_elpl = Lift / (2 * pi * rho * Uinf * le^2);

% ---- Creating the optimization equation ----
c = 2 * cos(fai) .* delta_sigma;
b = 3 * pi / 2 *(eta .* cos(fai) + zeta .* sin(fai)) .* delta_sigma;
for i = 1:N
	for j = 1:N
		A(i,j) = pi * q(i,j) .* delta_sigma(i);
	end
end

% ---- solve optimization problem ----
AAA = A + A';
ccc = -c;
bbb = -b;
LeftMat = vertcat([AAA ccc' bbb'], [ccc 0 0], [bbb 0 0]);
RightMat = vertcat(linspace(0,0,N)', -1, -beta);
SolveMat = LeftMat \ RightMat;
g = SolveMat(1:N);
mu = SolveMat(N+1:N+2);	% Lagrange multiplier


% ---- After Solve ----
% efficient : span efficiency
% Gamma : Local circulation
% InducedDrag : Sum of induced drag
% Vn : Local induced vertical velocisy
% Lift0_elpl : Lift at root culcurated by area of the ellipse circulation distribution
% Gamma0_elpl : Circulation at root
% Gamma_elpl : Analytical ellipse circulation distribution @ beta = 1
% ---------------------
efficiency_inverse = g' * A * g;
efficiency = 1 / efficiency_inverse;
Gamma = g * Lift / (2 * le * rho * Uinf);
InducedDrag = efficiency_inverse * InducedDrag_elpl;

Local_Lift = 4 * rho * Uinf * Gamma' .* cos(fai);
Vn = zeros(1,N);
for i = 1:N
	for j = 1:N
		Vn(i) =Vn(i)+ Q(i,j).*Gamma(j);
	end
end
Local_InducedDrag = rho * Gamma'.* Vn;

% ---- Aerodynamic force when Elliptical cierculation distribution----
Lift0_elpl = 2 * Lift / pi / le;
Lift_elpl = 4 * Lift0_elpl * sqrt(1 - (y / le).^2);
Gamma0_elpl = Lift0_elpl / (rho * Uinf * cos(fai(1)));
Gamma_elpl = Gamma0_elpl * sqrt(1 - (y / le).^2);
Local_InducedDrag_elpl = 2 * rho * Gamma_elpl * Vn_elpl;

% ---- Bending Moment ----
Local_BendingMoment = zeros(1,N);
Local_BendingMoment_elpl = zeros(1,N);
for i = 1:N
	tmp1 = 0;
	tmp2 = 0;
	for j = i:N
		tmp1 = tmp1 + Local_Lift(j) * (y(j) - y(i));
		tmp2 = tmp2 + Lift_elpl(j) * (y(j) - y(i));
	end
	Local_BendingMoment(i) = tmp1;
	Local_BendingMoment_elpl(i) = tmp2;
end

y = InducedDrag;

end