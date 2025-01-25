#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}


double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		int i = 0;
		double temp = 0;
		solution X0(x0), X1(x0 + d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if (X0.y == X1.y) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);

			return p;
		}
		
		if (X1.y > X0.y) {
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y)
			{
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x) - d;

				return p;
			}
		}
		do {
			if (solution::f_calls > Nmax) {
				X0.flag = 1;
				return 0;
			}	
			i = i + 1;
			temp = m2d(X0.x);
			X0 = X1;
			X1.x = x0 + pow(alpha, i) * d;
			X1.fit_fun(ff, ud1, ud2);
		} while (X0.y >= X1.y);
		if (d > 0) {
			p[0] = temp;
			p[1] = m2d(X1.x);
		}
		else {
			p[1] = temp;
			p[0] = m2d(X1.x);
		}

		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution a0(a), b0(b);
		a0.fit_fun(ff, ud1, ud2);
		b0.fit_fun(ff, ud1, ud2);
		solution c(0), d(0);
		int k = 1;
		while (GetFib(k) < (b0.x - a0.x) / epsilon) {
			k++;
		}
		c.x = b0.x - GetFib(k - 1) / GetFib(k) * (b0.x - a0.x);
		d.x = a0.x + b0.x - c.x;

		c.fit_fun(ff, ud1, ud2);
		d.fit_fun(ff, ud1, ud2);
		Xopt.ud = b - a;
		for (int i = 0; i < k - 3; i++) {

			if (c.y < d.y) {
				b0 = d;
			}
			else {
				a0 = c;
			}

			c.x = b0.x - GetFib(k - i - 2) / GetFib(k - i - 1) * (b0.x - a0.x);
			d.x = a0.x + b0.x - c.x;
			c.fit_fun(ff, ud1, ud2);
			d.fit_fun(ff, ud1, ud2);
			Xopt.ud.add_row(m2d(b0.x - a0.x));
		}
		Xopt = c;

		Xopt.flag = 0;

		return Xopt;
		
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int i = 0;
		double l, m;
		double c = (a + b) / 2.0;
		solution a0(a), b0(b), c0(c), d0(0), di(a);

		a0.fit_fun(ff, ud1, ud2);
		b0.fit_fun(ff, ud1, ud2);
		c0.fit_fun(ff, ud1, ud2);

		Xopt.ud = b - a;
		do {
			Xopt.ud.add_row(m2d(b0.x - a0.x));
			l = m2d(a0.y) * m2d(pow(b0.x, 2) - pow(c0.x, 2)) +
				m2d(b0.y) * m2d(pow(c0.x, 2) - pow(a0.x, 2)) +
				m2d(c0.y) * m2d(pow(a0.x, 2) - pow(b0.x, 2));

			m = m2d(a0.y) * m2d(b0.x - c0.x) +
				m2d(b0.y) * m2d(c0.x - a0.x) +
				m2d(c0.y) * m2d(a0.x - b0.x);

			if (m <= 0) {
				Xopt.flag = -1;
				return Xopt;
			}

			di.x = d0.x;
			d0.x = 0.5 * l / m;
			d0.fit_fun(ff, ud1, ud2);
			di.fit_fun(ff, ud1, ud2);

			if (a0.x < d0.x && d0.x < c0.x) {
				if (d0.y < c0.y) {
					b0 = c0;
					c0 = d0;
				}
				else {
					a0 = d0;
				}
			}
			else {
				if (c0.x < d0.x && d0.x < b0.x) {
					if (d0.y < c0.y)
					{
						a0 = c0;
						c0 = d0;
					}
					else {
						b0 = d0;
					}
				}
				else {
					Xopt = d0;
					Xopt.flag = -1;
					return Xopt;
				}
			}

			i = i + 1;
			Xopt.ud.add_row((b0.x - a0.x));

			if (solution::f_calls > Nmax) {
				Xopt.flag = 1;
				return 0;
			}

		} while ((b0.x - a0.x) >= epsilon && fabs(m2d(d0.x) - m2d(di.x)) >= gamma);

		Xopt = d0;
		Xopt.fit_fun(ff, ud1, ud2);

		Xopt.flag = 0;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}


solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution XB, XB_, X(x0);
		X.fit_fun(ff, ud1, ud2);
		matrix ud(x0);
		do
		{
			XB = X;
			X = HJ_trial(ff, XB, s, ud1, ud2);
			X.fit_fun(ff, ud1, ud2);
			if (X.y < XB.y) {
				do
				{
					XB_ = XB;
					XB = X;
					X.x = 2 * XB.x - XB_.x;
					X.fit_fun(ff, ud1, ud2);
					X = HJ_trial(ff, X, s, ud1, ud2);
					if (solution::f_calls > Nmax) {
						Xopt = XB;
						Xopt.flag = -1;
						return Xopt;
					}
				} while (X.y < XB.y);
				X = XB;
			}
			else
				s = alpha * s;
			if (solution::f_calls > Nmax) {
				Xopt = XB;
				Xopt.flag = -1;
				return Xopt;
			}
			ud.add_col(X.x);
		} while (s >= epsilon);
	
		Xopt = XB;
		Xopt.ud = ud;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_dim(XB);
		matrix e = ident_mat(n);
		solution X;
		for (int j = 0; j < n; j++) {
			X.x = XB.x + s * e[j];
			X.fit_fun(ff, ud1, ud2);

			if (X.y < XB.y) { 
				XB = X; 
			}
			else {
				X.x = XB.x - s * e[j];
				X.fit_fun(ff, ud1, ud2);

				if (X.y < XB.y) XB = X;
			}
		}

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		matrix dj, lambda, p, s, ud(x0);
		solution XB(x0), X;
		XB.fit_fun(ff, ud1, ud2);
		int n = get_dim(XB);
		dj = ident_mat(n);
		lambda = matrix(n, 1, 0.0);
		p = matrix(n, 1, 0.0);
		s = matrix(s0);
		int max_s;
		do
		{
			for (int j = 0; j < n; j++) {
				X.x = XB.x + s(j) * dj[j];
				X.fit_fun(ff, ud1, ud2);
				if (X.y < XB.y) {
					XB = X;
					lambda(j) = lambda(j) + s(j);
					s(j) = alpha * s(j);
				}
				else {
					s(j) = -beta * s(j);
					p(j) = p(j) + 1;
				}	
			}

			X = XB;

			bool reset = true;

			for (int j = 0; j < n; j++) {
				if (lambda(j) == 0 || p(j) == 0) {
					reset = false;
					break;
				}
			}

			if (reset) {
				matrix Q = dj;
				matrix v = matrix(n, n);
				matrix d_new = ident_mat(n);

				d_new = Q * d_new;
				v.set_col(d_new[0] / norm(d_new[0]), 0);
				for (int j = 0; j < n; j++) {
	
					matrix dot_product(n, 1);
					for (int k = 0; k < j; k++) {
						dot_product.set_col(dot_product[0] + (trans(d_new[j])) * dj[k] * dj[k], 0);
					}
					v.set_col((d_new[j] - dot_product[0]) / norm(d_new[j] - dot_product[0]), j);
					
				}

				dj = v;

				lambda = matrix(n, 1, 0.0);
				p = matrix(n, 1, 0.0);
				s = s0;
			}
			if (solution::f_calls > Nmax) {
				Xopt = X;
				Xopt.flag = -1;
				return Xopt;
			}
			ud.add_col(X.x);
			max_s = 0;
			for (int j = 1; j < n; j++) {
				if (abs(s(max_s)) < abs(s(j))) {
					max_s = j;
				}
			}
		} while (abs(s(max_s)) >= epsilon);

		Xopt = X;
		Xopt.ud = ud;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int min = 0;
		Xopt.flag = 0;

		int n = get_dim(x0);
	
		matrix p(n, n + 1);
		matrix e = ident_mat(n);

		p.set_col(x0, 0);
		for (int i = 1; i <= n; i++)
			p.set_col(p[0] + e[i - 1] * s, i); 

		solution::f_calls += 1 + n;
		matrix f(n + 1, 1);
		for (int i = 0; i <= n; i++)
			f(i) = m2d(ff(p[i], ud1, ud2));

		double max_diff;
		do {
			max_diff = 0.0;
			int max = 0;
			min = 0;
			for (int i = 1; i <= n; i++) {
				if (f(max) < f(i)) max = i;
				if (f(min) > f(i)) min = i;
			}
			if (max == min) {
				Xopt.flag = -1;
				break;
			}

			matrix p_(n, 1);
			for (int i = 0; i <= n; i++)
			{
				if (i != max) {
					p_.set_col(p_[0] + p[i], 0);
				}
			}

			p_.set_col(p_[0] / n, 0);
			matrix p_odb = p_[0] + (p_[0] - p[max]) * alpha;

			solution::f_calls++;
			double f_odb = m2d(ff(p_odb, ud1, ud2));

			if (f_odb < f(max))
			{
				matrix p_e = p_ + (p_odb[0] - p_[0]) * gamma;

				solution::f_calls++;
				double f_e = m2d(ff(p_e, ud1, ud2));

				if (f_e < f_odb)
				{
					p.set_col(p_e[0], max);
					f(max) = f_e;
				}
				else
				{
					p.set_col(p_odb[0], max);
					f(max) = f_odb;
				}
			}
			else
			{
				if (f(min) <= f_odb && f_odb < f(max))
				{
					p.set_col(p_odb[0], max);
					f(max) = f_odb;
				}
				else
				{
					matrix p_z = p_[0] + (p[max] - p_[0]) * beta;

					solution::f_calls++;
					double f_z = m2d(ff(p_z, ud1, ud2));

					if (f_z >= f(max))
					{
						for (int i = 0; i <= n; i++)
						{
							if (i != min) {
								p.set_col((p[i] + p[min]) * delta, i);
								solution::f_calls++;
								f(i) = m2d(ff(p[i], ud1, ud2));
							}
						}
					}
					else
					{
						p.set_col(p_z[0], max);
						f(max) = f_z;
					}
				}
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = -1;
				break;
			}
			for (int i = 0; i <= n; i++)
			{
				if (i != min) {
					double diff = norm(p[min] - p[i]);
					if (diff > max_diff)
						max_diff = diff;
				}
			}
		} while (max_diff > epsilon);
		Xopt.x = p[min];
		Xopt.y = f(min);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}


solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}


solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		int i = 0;

		double alpha = pow(static_cast<double>(N), -0.5);
		double beta = pow(2.0 * N, -0.25);

		// Zerujemy liczniki (o ile chcemy śledzić f_calls)
		solution::clear_calls();

		// ========================
		// Funkcje do losowania (lokalnie, bez makr)
		// ========================
		auto randU = []() -> double {
			// U(0,1)
			return (double)rand() / (double)RAND_MAX;
			};

		auto randN = [&]() -> double {
			// N(0,1) metodą Box-Muller (wersja uproszczona)
			double u1 = randU();
			double u2 = randU();
			return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.1415 * u2);
			};

		// ========================
		// 1) Inicjalizacja populacji
		// ========================
		vector<solution> P(mi);

		for (int j = 0; j < mi; j++) {
			// Tworzymy osobnika s z wektorem x (N×1)
			solution s(matrix(N, 1));
			// Rezerwujemy w s.ud miejsce na wektor sigma (N×1)
			s.ud = matrix(N, 1);

			// Losujemy każdą współrzędną x(d,0) w przedziale [ lb(d,0), ub(d,0) ]
			for (int d = 0; d < N; d++)
			{
				double r = randU();
				double xd = lb(d, 0) + r * (ub(d, 0) - lb(d, 0));
				s.x(d, 0) = xd;
				// Początkowa sigma
				s.ud(d, 0) = sigma0(d, 0);
			}

			// Wywołujemy funkcję celu, aby mieć y = f(x)
			s.fit_fun(ff, ud1, ud2);

			// Jeśli już na etapie inicjalizacji przekroczymy limit wywołań ff:
			if (solution::f_calls > Nmax) {
				throw string("EA: Przekroczono dozwoloną liczbę wywołań funkcji celu w trakcie inicjalizacji.");
			}
			P[j] = s;
		}

		// Pomocnicza lambda do znalezienia indeksu najlepszego (po s.y(0,0))
		auto bestIndex = [&](const vector<solution>& pop) {
			int bIdx = 0;
			double bVal = pop[0].y(0, 0);
			for (int j = 1; j < (int)pop.size(); j++)
			{
				double val = pop[j].y(0, 0);
				if (val < bVal)
				{
					bVal = val;
					bIdx = j;
				}
			}
			return bIdx;
			};

		// ========================
		// Główna pętla
		// ========================
		while (true)
		{
			// Sprawdzamy warunek stopu
			int bIdxPop = bestIndex(P);
			double fBest = P[bIdxPop].y(0, 0);

			// Czy f(x*) <= epsilon?
			if (fBest <= epsilon)
			{
				// Sukces
				P[bIdxPop].flag = 0;
				Xopt = P[bIdxPop];
				Xopt.flag = 0;
				return Xopt;
			}

			// Czy przekroczyliśmy Nmax wywołań funkcji celu?
			if (solution::f_calls >= Nmax)
			{
				// Koniec – limit wyczerpany
				Xopt = P[bIdxPop];
				Xopt.flag = -1;
				return Xopt;
			}

			// ========================
			// 2) Koło ruletki: fi_j = 1 / f_j
			// ========================
			vector<double> fi(mi), q(mi + 1, 0.0);
			double sumFi = 0.0;
			for (int j = 0; j < mi; j++)
			{
				double fj = P[j].y(0, 0);
				fi[j] = 1.0 / fj;
				sumFi += fi[j];
			}
			for (int j = 1; j <= mi; j++)
			{
				q[j] = q[j - 1] + (fi[j - 1] / sumFi);
			}

			// ========================
			// 3) Generowanie potomstwa (lambda)
			// ========================
			vector<solution> offspring(lambda);

			for (int off = 0; off < lambda; off++)
			{
				// Wybór rodzica A
				double rA = randU();
				int idxA = 0;
				for (int j = 1; j <= mi; j++)
				{
					if (rA <= q[j])
					{
						idxA = j - 1;
						break;
					}
				}
				// Wybór rodzica B
				double rB = randU();
				int idxB = 0;
				for (int j = 1; j <= mi; j++)
				{
					if (rB <= q[j])
					{
						idxB = j - 1;
						break;
					}
				}
				// Krzyżowanie: child.x(d,0) = A.x(d,0) + (1-rC)*B.x(d,0)
				double rC = randU();
				solution child(matrix(N, 1));
				child.ud = matrix(N, 1);

				for (int d = 0; d < N; d++)
				{
					double xA = P[idxA].x(d, 0);
					double xB = P[idxB].x(d, 0);
					double sA = P[idxA].ud(d, 0);
					double sB = P[idxB].ud(d, 0);

					child.x(d, 0) = xA + (1.0 - rC) * xB;
					// Sigma np. średnia
					child.ud(d, 0) = 0.5 * (sA + sB);
				}

				// Mutacja samoadaptacyjna
				double globalN = randN(); // globalne N(0,1)
				for (int d = 0; d < N; d++)
				{
					double localN = randN();
					double sigma_old = child.ud(d, 0);

					// sigma_d = sigma_d * exp(alpha*gN + beta*lN)
					double sigma_new = sigma_old * exp(alpha * globalN + beta * localN);
					child.ud(d, 0) = sigma_new;

					// x_d = x_d + sigma_d * N(0,1)
					double step = sigma_new * randN();
					child.x(d, 0) += step;

					// Granice
					if (child.x(d, 0) < lb(d, 0))
						child.x(d, 0) = lb(d, 0);
					if (child.x(d, 0) > ub(d, 0))
						child.x(d, 0) = ub(d, 0);
				}

				// Obliczamy f(child.x)
				child.fit_fun(ff, ud1, ud2);

				offspring[off] = child;
			}

			// ========================
		   // 4) Selekcja (μ+λ) -> μ
		   // ========================
			vector<solution> combined;
			combined.reserve(mi + lambda);

			// stare osobniki:
			for (int i = 0; i < mi; i++)
				combined.push_back(P[i]);
			// potomki
			for (int i = 0; i < lambda; i++)
				combined.push_back(offspring[i]);

			// Sortowanie rosnąco po y(0,0)
			std::sort(combined.begin(), combined.end(),
				[&](const solution& A, const solution& B) { return A.y(0, 0) < B.y(0, 0); });

			// Pierwsze mi -> nowa populacja
			for (int i = 0; i < mi; i++)
				P[i] = combined[i];

		} // koniec while

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}