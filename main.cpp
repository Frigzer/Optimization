/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab2();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	srand(time(NULL));

	double* ekspansja = new double[2];
	int Nmax = 1000;
	double epsilon = 1e-05, gamma = 1e-200;

	// Dane to tabeli 1
	ofstream exp_tab_1("./dane/lab_01/problem_testowy/exp_tab_1.txt");
	ofstream fib_tab_1("./dane/lab_01/problem_testowy/fib_tab_1.txt");
	ofstream lag_tab_1("./dane/lab_01/problem_testowy/lag_tab_1.txt");

	// Trzy wpó³czynniki alfa dla ekspansji
	double alpha, alpha_1 = 3.1, alpha_2 = 5.0, alpha_3 = 6.6;
	double x0;
	double d = 2.0;
	alpha = alpha_1;

	for (int i = 0; i < 300; i++)
	{
		x0 = -100 + (double)rand() / RAND_MAX * (200);

		if (i == 100) alpha = alpha_2;
		else if (i == 200) alpha = alpha_3;
		ekspansja = expansion(ff1T, x0, d, alpha, Nmax);
		exp_tab_1 << x0 << "\t" << ekspansja[0] << "\t" << ekspansja[1] << "\t" << solution::f_calls << endl;
		solution::clear_calls();

		solution fib1 = fib(ff1T, ekspansja[0], ekspansja[1], epsilon);
		fib_tab_1 << m2d(fib1.x) << "\t" << m2d(fib1.y) << "\t" << solution::f_calls << "\t" << fib1.flag << endl;
		solution::clear_calls();

		solution lag1 = lag(ff1T, ekspansja[0], ekspansja[1], epsilon, gamma, Nmax);
		lag_tab_1 << m2d(lag1.x) << "\t" << m2d(lag1.y) << "\t" << solution::f_calls << "\t" << lag1.flag << endl;
		solution::clear_calls();
	}
	exp_tab_1.close();
	fib_tab_1.close();
	lag_tab_1.close();

	// Dane do wykresu
	ofstream fib_wykres("./dane/lab_01/problem_testowy/fib_wykres.txt");
	ofstream lag_wykres("./dane/lab_01/problem_testowy/lag_wykres.txt");
	
	solution fib2 = fib(ff1T, -100, 100, epsilon);
	fib_wykres << fib2 << "\n\n" << fib2.ud << endl;
	solution::clear_calls();

	solution lag2 = lag(ff1T, -100, 100, epsilon, gamma, Nmax);
	lag_wykres << lag2 << "\n\n" << lag2.ud << endl;
	solution::clear_calls();

	fib_wykres.close();
	lag_wykres.close();

	// Dane do tabeli 3
	ofstream fib_tab_3("./dane/lab_01/problem_rzeczywisty/fib_tab_3.txt");
	ofstream lag_tab_3("./dane/lab_01/problem_rzeczywisty/lag_tab_3.txt");

	solution fib3 = fib(ff1R, 1e-4, 1e-2, epsilon);
	fib_tab_3 << fib3 << endl;
	solution::clear_calls();

	solution lag3 = lag(ff1R, 1e-4, 1e-2, epsilon, gamma, Nmax);
	lag_tab_3 << lag3 << endl;
	solution::clear_calls();

	fib_tab_3.close();
	lag_tab_3.close();

	// Dane do symulacji
	ofstream fib_sym("./dane/lab_01/problem_rzeczywisty/fib_symulacja.csv");
	ofstream lag_sym("./dane/lab_01/problem_rzeczywisty/lag_symulacja.csv");

	matrix Y0 = matrix(3, new double[3] {5, 1, 20});
	matrix* Y_fib = solve_ode(df1, 0, 1, 2000, Y0, NAN, fib3.x(0));
	solution::clear_calls();
	matrix* Y_lag = solve_ode(df1, 0, 1, 2000, Y0, NAN, lag3.x(0));

	fib_sym << Y_fib[1];
	lag_sym << Y_lag[1];

	fib_sym.close();
	lag_sym.close();
}

void lab2()
{
	double s = 0.77;
	double alphaHJ = 0.5;
	double aplhaR = 2;
	double beta = 0.5;
	double epsilon = 1e-3;
	int Nmax = 1000;
	solution opt;
	matrix x0, s0;
	s0 = matrix(2, 1, s);

	ofstream HookeToFile("hooke.txt");

	x0 = matrix(2, 1, 0.0);
	x0(0) = -0.3;
	x0(1) = 0;
	s = 0.1;
	s0 = matrix(2, 1, s);

	solution xopt = Rosen(ff2T, x0, s0, aplhaR, beta, epsilon, Nmax);
	cout << x0(0) << ";" << x0(1) << " " << xopt.x(0) << ";" << xopt.x(1) << " " << xopt.y;
	/*
	for (int i = 0; i < 100; i++)
	{
		x0 = 2 * rand_mat(2, 1) - 1;
		solution Hooke = HJ(ff2T, x0, s, alphaHJ, epsilon, Nmax);
		int a = solution::f_calls;
		HookeToFile << "x0: " << x0(0) << "x: " << Hooke.x << "y: " << Hooke.y << a << ";" << Hooke.flag << ";" << endl;
		x0(0);
		solution::clear_calls();
	}
	*/
}

void lab3()
{
	for (int i = 0; i < 100; i++)
	{

	}
}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
