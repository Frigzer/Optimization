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
		lab1();
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
	double epsilon = 1e-2, gamma = 1e-200;
	ofstream expToFile("./expansion.txt");
	ofstream fibToFile("./fibonacci.txt");
	ofstream lagToFile("./lagrange.txt");

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
		expToFile << x0 << "\t" << ekspansja[0] << "\t" << ekspansja[1] << "\t" << solution::f_calls << endl;
		solution::clear_calls();

		solution fibonacci = fib(ff1T, ekspansja[0], ekspansja[1], epsilon);
		fibToFile << m2d(fibonacci.x) << "\t" << m2d(fibonacci.y) << "\t" << solution::f_calls << "\t" << fibonacci.flag << endl;
		solution::clear_calls();

		solution lagrange = lag(ff1T, ekspansja[0], ekspansja[1], epsilon, gamma, Nmax);
		lagToFile << m2d(lagrange.x) << "\t" << m2d(lagrange.y) << "\t" << solution::f_calls << "\t" << lagrange.flag << endl;
		solution::clear_calls();
		x0++;
	}

	ofstream fibWykres("./fibWykres.txt");
	ofstream lagWykres("./lagWykres.txt");
	
	solution opt = fib(ff1T, -100, 100, epsilon);
	fibWykres << opt << "\n\n" << opt.ud << endl;

	solution opt2 = lag(ff1T, -100, 100, epsilon, gamma, Nmax);
	lagWykres << opt2 << "\n\n" << opt2.ud << endl;


}

void lab2()
{

}

void lab3()
{

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
