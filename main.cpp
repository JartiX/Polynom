#include <iostream>
#include "polynom/polynom.h"
#include "Rational/Rational.h"

using namespace std;

int main() {
    Polynomial<double> pol({ 3,3,4,9 });
    Polynomial<double> poly({ 2,3,5, 8, 6});
    cout << pol << "POWER " << pol.GetPower()<<endl;
    cout << poly << endl;

    cout << pol * poly;
    cout << "pol / 2: " << poly / 2 << endl;

    cout << "pol / poly: " << pol / poly << endl;

    cout << pol.integrate() << endl;
    cout << poly.derivative() << endl;

    Polynomial<Rational> rat({ Rational(3,2), Rational(2,5), Rational(6,4) });

    cout << rat << endl;

    cout << rat.derivative() << endl;
    cout << rat.integrate() << endl;
    cout << rat.evaluate(Rational(2,8)) << endl;
    cout << Rational(2, 8) * Rational(2, 8) << endl;
    cout << Rational(2, 8) + Rational(3, 2);

    /*Polynomial<double> poll;
    cin >> poll;
    cout << poll;*/

}