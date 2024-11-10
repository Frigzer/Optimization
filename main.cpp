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
	srand(time(NULL));

	double alpha = 2;
	double epsilon = 1e-06;
	int Nmax = 1000;

	// Dane do tabeli 1
	ofstream x_val_tab_1("./dane/lab_02/problem_testowy/x_val_tab_1.txt");
	ofstream hooke_tab_1("./dane/lab_02/problem_testowy/hooke_tab_1.txt");
	ofstream rosen_tab_1("./dane/lab_02/problem_testowy/rosen_tab_1.txt");

	// Trzy dlugosci kroku
	double s, s_1 = 0.7, s_2 = 0.45, s_3 = 0.1;
	double beta = 0.5;
	s = s_1;
	matrix x(2, 1);
	
	bool wykres = false;

	solution hooke2, rosen2;
	for (int i = 0; i < 300; i++)
	{
		x(0) = -1 + static_cast<double>(rand()) / RAND_MAX * (2);
		x(1) = -1 + static_cast<double>(rand()) / RAND_MAX * (2);

		if (i == 100) s = s_2;
		else if (i == 200) s = s_3;

		x_val_tab_1 << x(0) << "\t" << x(1) << endl;

		solution hooke1 = HJ(ff2T, x, s, beta, epsilon, Nmax);
		hooke_tab_1 << m2d(hooke1.x(0)) << "\t" << m2d(hooke1.x(1)) << "\t" << m2d(hooke1.y) << "\t" << solution::f_calls << "\t" << endl;
		solution::clear_calls();

		solution rosen1 = Rosen(ff2T, x, matrix(2, 1, s), alpha, beta, epsilon, Nmax);
		rosen_tab_1 << m2d(rosen1.x(0)) << "\t" << m2d(rosen1.x(1)) << "\t" << m2d(rosen1.y) << "\t" << solution::f_calls << "\t" << endl;
		solution::clear_calls();

		if (!wykres && (abs(m2d(hooke1.y)) < epsilon) && (abs(m2d(rosen1.y)) < epsilon)) {
			wykres = true;
			hooke2 = hooke1;
			rosen2 = rosen1;
		}
	}
	x_val_tab_1.close();
	hooke_tab_1.close();
	rosen_tab_1.close();
	
	// Dane do wykresu
	ofstream hooke_wykres("./dane/lab_02/problem_testowy/hooke_wykres.txt");
	ofstream rosen_wykres("./dane/lab_02/problem_testowy/rosen_wykres.txt");


	int* hooke2_size = get_size(hooke2.ud);
	for (int i = 0; i < hooke2_size[1]; i++) {
		hooke_wykres << hooke2.ud(0, i) << "\t" << hooke2.ud(1, i) << endl;
	}
	solution::clear_calls();

	
	int* rosen2_size = get_size(rosen2.ud);
	for (int i = 0; i < rosen2_size[1]; i++) {
		rosen_wykres << rosen2.ud(0, i) << "\t" << rosen2.ud(1, i) << endl;
	}
	solution::clear_calls();
	
	delete[] hooke2_size;
	delete[] rosen2_size;
	
	hooke_wykres.close();
	rosen_wykres.close();

	// Dane do tabeli 3
	ofstream hooke_tab_3("./dane/lab_02/problem_rzeczywisty/hooke_tab_3.txt");
	ofstream rosen_tab_3("./dane/lab_02/problem_rzeczywisty/rosen_tab_3.txt");

	s = 0.5;
	double k[2] = { 4.0, 8.0 };
	matrix x0 = matrix(2, k);

	solution hooke3 = HJ(ff2R, x0, s, beta, epsilon, Nmax);
	hooke_tab_3 << m2d(hooke3.x(0)) << "\t" << m2d(hooke3.x(1)) << "\t" << m2d(hooke3.y) << "\t" << solution::f_calls << "\t" << hooke3.flag << endl;
	solution::clear_calls();

	solution rosen3 = Rosen(ff2R, x0, matrix(2, 1, s), alpha, beta, epsilon, Nmax);
	rosen_tab_3 << m2d(rosen3.x(0)) << "\t" << m2d(rosen3.x(1)) << "\t" << m2d(rosen3.y) << "\t" << solution::f_calls << "\t" << rosen3.flag << endl;
	solution::clear_calls();

	hooke_tab_3.close();
	rosen_tab_3.close();

	// Dane do symulacji
	ofstream hooke_sym("./dane/lab_02/problem_rzeczywisty/hooke_symulacja.csv");
	ofstream rosen_sym("./dane/lab_02/problem_rzeczywisty/rosen_symulacja.csv");

	matrix Y0 = matrix(2, new double[2] {0.0, 0.0});

	matrix* Y_hooke = solve_ode(df2, 0, 0.1, 100, Y0, hooke3.x(0), hooke3.x(1));
	solution::clear_calls();

	matrix* Y_rosen = solve_ode(df2, 0, 0.1, 100, Y0, rosen3.x(0), rosen3.x(1));
	solution::clear_calls();

	hooke_sym << Y_hooke[1];
	rosen_sym << Y_rosen[1];

	hooke_sym.close();
	rosen_sym.close();
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
