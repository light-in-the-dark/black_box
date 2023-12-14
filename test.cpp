#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <sstream>
#include <cctype>

double eps = 1e-6;

template <typename T, typename P> double neg_sqrt(T val, P power = 1.0) {
    return val >= 0 ? pow(val, P(1.) / power) : -pow(-val, P(1.) / power);
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class Test {
private:
    std::stringbuf sbuf;
    std::streambuf* oldbuf;
    int n;

    void update_state() {
        this->n += 1;
        this->sbuf = std::stringbuf(std::ios::out);
        this->oldbuf = std::cout.rdbuf(std::addressof(this->sbuf));
    }
public:
    Test() {
        this->update_state();
    };
    void out_stream_to_be(std::string expect) {

        std::cout.rdbuf(this->oldbuf);

        std::string output = this->sbuf.str();

        if (output == expect) {
            std::cout << "Тест " << this->n << " пройден\n";
        }
        else {
            std::cout << "Тест " << this->n << " не пройден.\n" << "Ожидаемый результат: " << expect << "\n" << "Полученный результат: " << output << "\n";
        }
        this->update_state();
    }
};

int solve_linear(double a, double b) {
    double x = -b / a;
    std::cout << x << std::endl;
    return 0;
}

int solve_quadratic_equation(double a, double b, double c) {
    double D = std::pow(b, 2) - 4 * a * c;
    if (D < 0) {
        std::cout << "Нет вещественных корней" << std::endl;
        return 0;
    }
    if (D == 0) {
        double x = -b / (2 * a);
        std::cout << x << std::endl;
        return 0;
    }
    double x1 = (-b + sqrt(D)) / (2 * a);
    double x2 = (-b - sqrt(D)) / (2 * a);
    std::cout << x1 << " " << x2 << std::endl;
    return 0;
}

int solve_cubic_equation(double del, double a, double b, double c)
{
    a /= del;
    b /= del;
    c /= del;
    del /= del;

    double Q = (pow(a, 2.) - 3. * b) / 9.;
    double R = (2 * pow(a, 3.) - 9. * a * b + 27. * c) / 54.;

    double S = pow(Q, 3.) - pow(R, 2.);

    if (S > 0) {
        double sqrt_Q = neg_sqrt(Q, 2.);
        double phi = 1. / 3. * acos(R / sqrt(pow(Q, 3.)));
        double x1 = -2 * sqrt_Q * cos(phi) - a / 3;
        double x2 = -2 * sqrt_Q * cos(phi + 2. / 3. * M_PI) - a / 3.;
        double x3 = -2 * sqrt_Q * cos(phi - 2. / 3. * M_PI) - a / 3.;
        std::cout << x1 << " " << x2 << " " << x3 << std::endl;
        return 0;
    }

    if (std::abs(S) <= eps) {
        double sqrt_quibic_R = neg_sqrt(R, 3.);
        double x1 = -2. * sqrt_quibic_R - a / 3.;
        double x2 = sqrt_quibic_R - a / 3.;
        std::cout << x1 << " " << x2 << std::endl;
        return 0;
    }
    if (Q > 0) {
        double sqrt_Q = neg_sqrt(Q, 2.);
        double sqrt_quibic_Q = neg_sqrt(pow(Q, 3), 2.);
        double phi = 1. / 3. * acosh(abs(R) / sqrt_quibic_Q);
        double x = -2. * sgn(R) * sqrt_Q * cosh(phi) - a / 3.;
        std::cout << x << std::endl;
        return 0;
    }
    if (Q < 0) {
        double sqrt_quibic_Q = neg_sqrt(pow(abs(Q), 3), 2.);
        double phi = 1. / 3. * asinh(abs(R) / sqrt_quibic_Q);
        double x = -2. * sgn(R) * sqrt(abs(Q)) * sinh(phi) - a / 3.;
        std::cout << x << std::endl;
        return 0;
    }
    double x = -neg_sqrt(c - pow(a, 3.) / 27., 3.) - a / 3.;
    std::cout << x << std::endl;
    return 0;
}

double solve_equation(double a, double b, double c, double d) {
    if (a != 0) {
        return solve_cubic_equation(a, b, c, d);
    }
    if (b != 0) {
        return solve_quadratic_equation(b, c, d);
    }
    if (c != 0) {
        return solve_linear(c, d);
    }
    if (d != 0) {
        std::cout << "Нет решений" << std::endl;
        return 0;
    }

    std::cout << "Любое значение" << std::endl;

    return 0;
}

int main() {
    Test test = Test();
    setlocale(LC_ALL, "Russian");
    solve_equation(1, -3, -6, 8);
    test.out_stream_to_be("-2 4 1\n");
    solve_equation(1, -9, 24, -16);
    test.out_stream_to_be("1 4\n");
    solve_equation(2, -11, 12, 19);
    test.out_stream_to_be("-0.839222\n");
    solve_equation(0, 1, -5, 4);
    test.out_stream_to_be("4 1\n");
    solve_equation(0, 9, -6, 1);
    test.out_stream_to_be("0.333333\n");
    solve_equation(0, 9, 1, 9);
    test.out_stream_to_be("Нет вещественных корней\n");
    solve_equation(0, 0, 1, 9);
    test.out_stream_to_be("-9\n");
    solve_equation(0, 0, 0, 9);
    test.out_stream_to_be("Нет решений\n");
    solve_equation(0, 0, 0, 0);
    test.out_stream_to_be("Любое значение\n");
}

