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
		//int i = 0;
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
		//Tu wpisz kod funkcji

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
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
