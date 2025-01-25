#include"user_funs.h"

const double PI = 3.1415; // Wybranie tego co dziala



matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	double PI = 3.14;
	y = -cos(0.1 * m2d(x)) * exp(-pow(0.1 * m2d(x) - 2 * PI, 2)) + 0.002 * pow(0.1 * x, 2);
	return y;
}

double GetFib(int n) {
	if (n <= 1) return n;
	double a = 0, b = 1, c;
	for (int i = 2; i <= n; ++i) {
		c = a + b;
		a = b;
		b = c;
	}
	return b;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(3, 1);
	double a = 0.98, b = 0.63, g = 9.81;

	double PA = 0.5;
	double TA = 90.0;

	double PB = 1.0;
	double DB = 0.00365665;

	double Fin = 0.01;
	double Tin = 20.0;

	double Aout = Y(0) > 0 ? a * b * m2d(ud2) * sqrt(2 * g * Y(0) / PA) : 0;
	double Bout = Y(1) > 1 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;

	dY(0) = -Aout;
	dY(1) = Aout + Fin - Bout;
	dY(2) = Fin / Y(1) * (Tin - Y(2)) + Aout / Y(1) * (TA - Y(2));

	return dY;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});

	double t0 = 0.0;
	double tend = 2000.0;
	double dt = 1.0;

	matrix* Y = solve_ode(df1, t0, dt, tend, Y0, ud1, x);
	double max = Y[1](0, 2);

	int length = get_len(Y[0]);
	double y;
	for (int i = 1; i < length; i++)
	{
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	}
	y = abs(max - 50);
	return y;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	return matrix(pow(x(0)) + pow(x(1)) - cos(2.5 * PI * x(0)) - cos(2.5 * PI * x(1)) + 2);
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);

	double l = 1.0;
	double mr = 1.0;
	double mc = 5.0;
	double b = 0.5;
	double I = (mr / 3 + mc) * pow(l, 2);

	double k1 = m2d(ud1);
	double k2 = m2d(ud2);

	double alpha_ref = PI;
	double omega_ref = 0.0;

	double alpha = alpha_ref - Y(0);
	double omega = omega_ref - Y(1);
	double Mt = k1 * alpha + k2 * omega;

	dY(0) = Y(1);
	dY(1) = (Mt - b * Y(1)) / I;

	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	matrix Y0 = matrix(2, new double[2]{0, 0});

	matrix k1 = x(0);
	matrix k2 = x(1);

	double t0 = 0.0;
	double tend = 100.0;
	double dt = 0.1;

	matrix* Y = solve_ode(df2, t0, dt, tend, Y0, k1, k2);
	int length = get_len(Y[0]);

	double Q = 0.0;
	for (int i = 0; i < length; i++)
	{
		double alpha = Y[1](i, 0) - PI;
		double omega = Y[1](i, 1);

		double Mt = m2d(k1) * alpha + m2d(k2) * omega;

		Q += (10 * pow(alpha, 2) + pow(omega, 2) + pow(Mt, 2)) * dt;
	}

	return matrix(Q);
}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
	double y = PI * std::sqrt(m2d(pow(x(0) * PI, 2) + pow(x(1) * PI, 2)));
	return sin(y) / y;

}

matrix ff4T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	cout << "breaks1";
	x = ud2;
	if (isnan(ud2(0, 0)))
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	else
		y = ff4T(ud2[0] + x * ud2[1], 0, ud1);
	return y;
}

matrix gf4T(matrix x, matrix ud1, matrix ud2) {
	matrix g(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
	return g;
}
matrix hf4T(matrix x, matrix ud1, matrix ud2) {
	matrix H(2, 2);
	H(0, 0) = 10;
	H(0, 1) = 8;
	H(1, 0) = 8;
	H(1, 1) = 10;
	return H;
}
}

matrix ff6T(matrix x, matrix ud1, matrix ud2) {
	static double d_2o5_pi = 2.5 * _Pi_val;
	return x(0) * x(0) + x(1) * x(1) - cos(d_2o5_pi * x(0)) - cos(d_2o5_pi * x(1)) + 2.0;
}

matrix df6(double t, matrix Y, matrix ud1, matrix ud2)
{
	double m1 = 5.0, m2 = 5.0;
	double k1 = 1.0, k2 = 1.0;
	double F = 1.0;

	double b1 = ud1(0, 0);
	double b2 = ud1(1, 0);

	double x1 = Y(0, 0);
	double v1 = Y(1, 0);
	double x2 = Y(2, 0);
	double v2 = Y(3, 0);

	double a1 = (-b1 * v1 - b2 * (v1 - v2) - k1 * x1 - k2 * (x1 - x2)) / m1;
	double a2 = (F + b2 * (v1 - v2) + k2 * (x1 - x2)) / m2;

	matrix dY(4, 1);
	dY(0, 0) = v1; // dx1/dt = v1
	dY(1, 0) = a1; // dv1/dt = a1
	dY(2, 0) = v2; // dx2/dt = v2
	dY(3, 0) = a2; // dv2/dt = a2

	return dY;
}

matrix ff6R(matrix x, matrix ud1, matrix ud2)
{
	double b1 = x(0, 0);
	double b2 = x(1, 0);
	matrix experimentalData = ud2;

	double dt = 0.1;
	double T = 100.0;
	double t0 = 0.0;

	// Pocz¹tkowe wartoœci
	matrix Y0(4, 1);
	Y0(0, 0) = 0.0; // pozycja x1
	Y0(1, 0) = 0.0; // prêdkoœæ v1
	Y0(2, 0) = 0.0; // pozycja x2
	Y0(3, 0) = 0.0; // prêdkoœæ v2

	matrix* S = solve_ode(df6, t0, dt, T, Y0, x);

	// Wyniki równania ró¿niczkowego
	matrix time = S[0];      // Kroki czasowe
	matrix positions = S[1]; // Pozycje x1, x2

	int numPoints = get_len(time);

	// Liczenie b³êdu dopasowania
	double error = 0.0;
	for (int i = 0; i < numPoints; ++i)
	{
		double sim_x1 = positions(i, 0);
		double sim_x2 = positions(i, 2);
		double exp_x1 = experimentalData(i, 0);
		double exp_x2 = experimentalData(i, 1);

		error += pow(sim_x1 - exp_x1, 2) + pow(sim_x2 - exp_x2, 2);
	}

	// Sprawdzamy, czy aktualna para b1, b2 jest najlepsza
	static double best_error = std::numeric_limits<double>::max();
	static matrix best_positions = positions;
	static matrix best_time = time;
	static double best_b1 = 0, best_b2 = 0;

	if (error < best_error)
	{
		best_error = error;
		best_positions = positions;
		best_time = time;
		best_b1 = b1;
		best_b2 = b2;
	}

	std::ofstream file("output/lab6/wyniki_symulacji_p6.csv");
	if (!file.is_open())
	{
		std::cerr << "B³¹d otwierania pliku do zapisu!" << std::endl;
		return matrix(1, 1, 1e9); // B³¹d zwracany jako du¿a wartoœæ
	}

	file << "t[s],x1,x2\n";
	for (int i = 0; i < get_len(best_time); ++i)
	{
		file << best_time(i, 0) << "," << best_positions(i, 0) << "," << best_positions(i, 2) << "\n";
	}

	file.close();

	delete[] S;

	matrix result(1, 1);
	result(0, 0) = best_error;
	return result;
}