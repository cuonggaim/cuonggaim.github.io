
#include <vector>
#include <cstdarg>
#include <climits>
#include <iostream>
#include <cmath>

#define T_INITIAL 1.0
#define ALPHA 0.9
#define T_MIN 0.00001
#define ITERATIONS 200

using namespace std;

class Point
{
private:
	vector<double> values;
public:
	Point(vector<double> vals)
	{
		values = vals;
	}

	vector<double> get_values() const
	{
		return values;
	}
};

class Domain
{
private:
	double lower_bound;
	double upper_bound;
public:
	Domain(double low, double high)
	{
		lower_bound = low;
		upper_bound = high;
	}

	double get_lower_bound() const
	{
		return lower_bound;
	}

	double get_upper_bound() const
	{
		return upper_bound;
	}
};

class MathFunction
{
private:
	double(*function)(Point&);
	vector<Domain> domains;
	int variables;
public:

	void set_function(double(*func)(Point&))
	{
		function = func;
	}

	void set_variables(int vars)
	{
		variables = vars;
	}

	void set_domains(vector<Domain> d)
	{
		domains = d;
	}

	double operator()(double x, ...)
	{

		vector<double> vals;
		vals.push_back(x);
		va_list ap;
		va_start(ap, x);
		for (int i = 0; i < variables - 1; ++i)
		{
			double tmp = va_arg(ap, double);
			vals.push_back(tmp);
		}
		va_end(ap);

		return (*function)(Point(vals));
	}

	double operator()(Point& point)
	{
		return (*function)(point);
	}

	vector<Domain> get_domains() const
	{
		return domains;
	}

	int get_variables() const
	{
		return variables;
	}
};

static float generate_rand_double(double low, double high)
{
	return low + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (high - low)));
}

Point get_initial_solution(MathFunction& f)
{
	vector<double> vals;
	for (int i = 0; i < f.get_variables(); ++i)
	{
		double lower_bound = (f.get_domains()[i]).get_lower_bound();
		double upper_bound = (f.get_domains()[i]).get_upper_bound();

		vals.push_back(generate_rand_double(lower_bound, upper_bound));
	}

	return Point(vals);
}

bool accept(double old_cost, double new_cost, double temperature)
{
	return exp((old_cost - new_cost) / temperature);
}

Point simmulated_annealing(MathFunction& f)
{
	Point solution = get_initial_solution(f);
	double old_cost = f(solution);
	double temperature = T_INITIAL;

	while (temperature > T_MIN)
	{
		for (int i = 0; i < ITERATIONS; ++i)
		{
			Point new_solution = get_initial_solution(f);
			double new_cost = f(new_solution);

			if (accept(old_cost, new_cost, temperature) > rand() % 2)
			{
				solution = new_solution;
				old_cost = new_cost;
			}
		}
		temperature *= ALPHA;
	}

	return solution;
}


template<int n>
double beale_func(Point& p)
{
	double x = p.get_values()[0];
	double y = p.get_values()[1];

	return pow((1.5 - x + x * y), 2.0) +
		(2.25 - x + x * pow(y, 2.0)) +
		pow((2.625 - x + x * pow(y, 3.0)), 2.0);
}

template<int n>
double sum_squares_func(Point& p)
{
	double result = 0;
	for (int i = 1; i <= n; ++i)
	{
		result += i * pow(p.get_values()[i - 1], 2.0);
	}

	return result;
}

int main()
{
	// Computing the minimum for the beale 
	MathFunction beale;
	beale.set_function(beale_func<2>);
	beale.set_variables(2);
	vector<Domain> dm;
	dm.push_back(Domain(-4.5, 4.5));
	dm.push_back(Domain(-4.5, 4.5));
	beale.set_domains(dm);

	Point average_beale(vector<double> {0.0, 0.0});
	Point max_beale(vector<double> {-DBL_MAX, -DBL_MAX});
	Point min_beale(vector<double> {DBL_MAX, DBL_MAX});

	for (int i = 0; i < 30; ++i)
	{
		Point p = simmulated_annealing(beale);
		vector<double> g = p.get_values();
		vector<double> current_sol = average_beale.get_values();

		average_beale = Point(vector<double> {g[0] + current_sol[0], g[1] + current_sol[1]});

		if (g[0] >= max_beale.get_values()[0] && g[1] >= max_beale.get_values()[1])
		{
			max_beale = Point(vector<double> {g[0], g[1]});
		}

		if (g[0] <= min_beale.get_values()[0] && g[1] <= min_beale.get_values()[1])
		{
			min_beale = Point(vector<double> {g[0], g[1]});
		}
	}


	cout << "The average solution is "
		<< (average_beale.get_values()[0] / 200) << " "
		<< (average_beale.get_values()[1] / 200)
		<< endl;

	cout << "The minimum solution found is"
		<< (min_beale.get_values()[0]) << " "
		<< (min_beale.get_values()[1])
		<< endl;

	cout << "The maximum solution found is "
		<< (max_beale.get_values()[0]) << " "
		<< (max_beale.get_values()[1])
		<< endl;

	// Computing the mnimum for the sum squares function
	MathFunction sum_squares;
	sum_squares.set_function(sum_squares_func<3>);
	sum_squares.set_variables(3);
	vector<Domain> dm2;
	dm2.push_back(Domain(-10, 10));
	dm2.push_back(Domain(-10, 10));
	dm2.push_back(Domain(-10, 10));
	sum_squares.set_domains(dm2);

	Point average_sol(vector<double> {0.0, 0.0, 0.0});
	Point max_sol(vector<double> {-DBL_MAX, -DBL_MAX, -DBL_MAX});
	Point min_sol(vector<double> {DBL_MAX, DBL_MAX, DBL_MAX});

	for (int i = 0; i < 30; ++i)
	{
		Point& p2 = simmulated_annealing(sum_squares);

		vector<double> g2 = p2.get_values();
		vector<double> current_sol = average_sol.get_values();

		average_sol = Point(vector<double> {g2[0] + current_sol[0],
			g2[1] + current_sol[1], g2[2] + current_sol[2]});

		if (g2[0] >= max_sol.get_values()[0]
			&& g2[1] >= max_sol.get_values()[1]
			&& g2[2] >= max_sol.get_values()[2])
		{
			max_sol = Point(vector<double> {g2[0], g2[1], g2[2]});
		}

		if (g2[0] <= min_sol.get_values()[0]
			&& g2[1] <= min_sol.get_values()[1]
			&& g2[2] <= min_sol.get_values()[2])
		{
			max_sol = Point(vector<double> {g2[0], g2[1], g2[2]});
		}
	}

	cout << "The average solution is: "
		<< (average_sol.get_values()[0] / 200) << " "
		<< (average_sol.get_values()[1] / 200) << " "
		<< (average_sol.get_values()[2] / 200) << " "
		<< endl;

	cout << "The minimum sollution found is "
		<< (min_sol.get_values()[0]) << " "
		<< (min_sol.get_values()[1]) << " "
		<< (min_sol.get_values()[2]) << " "
		<< endl;

	cout << "The maximum solution found so far is "
		<< (max_sol.get_values()[0]) << " "
		<< (max_sol.get_values()[1]) << " "
		<< (max_sol.get_values()[2]) << " "
		<< endl;

    return 0;
}
