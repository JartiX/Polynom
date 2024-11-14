#pragma once
#ifndef POLYNOM_H
#define POLYNOM_H
#include <iostream>
#include <vector>

using namespace std;

template <class T>
class Polynomial {
private:
    size_t power;             // Степень многочлена
    vector<T> coefficients;   // Коэффициенты многочлена
    void update_coef();

public:
    // Конструктор
    Polynomial();
    Polynomial(const vector<T>& coef);

    size_t GetPower();


    // Функция для нахождения значения многочлена при заданном x
    T evaluate(T x) const;

    Polynomial& operator =(const Polynomial<T>& r);
    Polynomial& operator +=(const Polynomial<T>& r);
    Polynomial& operator -=(const Polynomial<T>& r);
    Polynomial& operator *=(const Polynomial<T>& r);
    Polynomial& operator /=(const Polynomial& divisor);
    Polynomial& operator /=(const T& r);

    Polynomial operator +(const Polynomial<T>& r) const;
    Polynomial operator -(const Polynomial<T>& r) const;
    Polynomial operator *(const Polynomial<T>& r) const;
    Polynomial operator /(const Polynomial& divisor) const;
    Polynomial operator /(const T& r) const;

    Polynomial integrate();
    Polynomial derivative(int order = 1);

    // Вывод многочлена
    template <class T>
    friend std::ostream& operator<<(std::ostream& os, const Polynomial<T>& poly);

    // Ввод многочлена
    template <class T>
    friend std::istream& operator>>(std::istream& in, Polynomial<T>& r);
};


template <class T>
Polynomial<T>::Polynomial() {
    power = 0;
}

// Конструктор
template <class T>
Polynomial<T>::Polynomial(const vector<T>& coef) {
    coefficients = coef;
    update_coef();
}

template <class T>
size_t Polynomial<T>::GetPower() {
    return power;
}
template <class T>
void Polynomial<T>::update_coef() {
    power = coefficients.size() - 1;
}

// Функция для нахождения значения многочлена при заданном x
template <class T>
T Polynomial<T>::evaluate(T x) const {
    T result = 0;
    T powered_x = 1;
    for (int i = power; i >= 0; i--) {
        if (i != power) {
            powered_x *= x;
        }
        result = result + powered_x * coefficients[i];
    }
    return result;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator =(const Polynomial<T>& r) {
    coefficients = r.coefficients;
    power = r.power;
    return *this;
}

// Сложение двух многочленов
template <class T>
Polynomial<T> Polynomial<T>::operator +(const Polynomial<T>& other) const {
    Polynomial<T> copy_poly(*this);
    return copy_poly += other;

}

template <class T>
Polynomial<T>& Polynomial<T>::operator +=(const Polynomial<T>& other) {
    int maxDegree = max(power, other.power);
    vector<T> resultCoeffs(maxDegree + 1, 0);

    for (int i = 0; i <= power; ++i) {
        resultCoeffs[i] += coefficients[i];
    }
    for (int i = 0; i <= other.power; ++i) {
        resultCoeffs[i] += other.coefficients[i];
    }
    coefficients = resultCoeffs;
    update_coef();
    return *this;
}

// Вычитание двух многочленов
template <class T>
Polynomial<T> Polynomial<T>::operator -(const Polynomial<T>& other) const {
    Polynomial<T> copy_poly(*this);
    return copy_poly -= other;

}

template <class T>
Polynomial<T>& Polynomial<T>::operator -=(const Polynomial<T>& other) {
    int maxDegree = max(power, other.power);
    vector<T> resultCoeffs(maxDegree + 1, 0);

    for (int i = 0; i <= power; ++i) {
        resultCoeffs[i] += coefficients[i];
    }
    for (int i = 0; i <= other.power; ++i) {
        resultCoeffs[i] -= other.coefficients[i];
    }
    coefficients = resultCoeffs;
    update_coef();
    return *this;
}


// Умножение двух многочленов
template <class T>
Polynomial<T> Polynomial<T>::operator *(const Polynomial<T>& other) const {
    Polynomial<T> copy_pol(*this);
    return copy_pol *= other;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator *=(const Polynomial<T>& other) {
    int resultSize = power + other.power;
    vector<T> result(resultSize + 1, 0);

    for (int i = 0; i <= power; ++i) {
        for (int j = 0; j <= other.power; ++j) {
            result[i + j] += coefficients[i] * other.coefficients[j];
        }
    }
    coefficients = result;
    update_coef();

    return *this;
}

// Деление многочленов (деление на скаляр)
template <class T>
Polynomial<T> Polynomial<T>::operator /(const T& r) const {
    Polynomial<T> copy_pol(*this);
    return copy_pol /= r;
}
template <class T>
Polynomial<T>& Polynomial<T>::operator /=(const T& r) {
    for (int i = 0; i < coefficients.size(); ++i) {
        coefficients[i] /= r;
    }
    return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator /=(const Polynomial<T>& divisor) {
    if (divisor.power == 0 && divisor.coefficients[0] == 0) {
        throw runtime_error("Division by zero polynomial");
    }

    if (power < divisor.power) {
        coefficients = vector<T>(1, 0);
        power = 0;
        return *this;
    }

    vector<double> quotient_coeffs(power - divisor.power + 1, 0);
    vector<T> remainder = coefficients;

    while (remainder.size() >= divisor.coefficients.size()) {
        // Вычисляем текущий член частного
        size_t degree_diff = remainder.size() - divisor.coefficients.size();
        T coeff = remainder.back() / divisor.coefficients.back();
        quotient_coeffs[degree_diff] = coeff;

        // Вычитаем произведение текущего члена частного на делитель
        for (int i = 0; i < divisor.coefficients.size(); ++i) {
            remainder[i + degree_diff] -= coeff * divisor.coefficients[i];
        }

        // Удаляем ведущие нули из остатка
        while (!remainder.empty() && abs(remainder.back()) < 1e-10) {
            remainder.pop_back();
        }
    }

    coefficients = quotient_coeffs;
    update_coef();
    return *this;
}

template <class T>
Polynomial<T> Polynomial<T>::operator /(const Polynomial<T>& divisor) const {
    Polynomial<T> result(*this);
    return result /= divisor;
}

//производная
template <class T>
Polynomial<T> Polynomial<T>::derivative(int order) {
    vector<T> result = coefficients;

    for (int i = 0; i < order; ++i) {
        if (result.size() <= 1) {
            result.clear();
            break;
        }

        vector<T> der_coeffs(result.size() - 1);
        for (int j = 1; j < result.size(); ++j) {
            der_coeffs[j - 1] = result[j] * j;
        }
        result = der_coeffs;
    }
    coefficients = result;
    update_coef();
    return *this;
}

// Интеграл многочлена
template <class T>
Polynomial<T> Polynomial<T>::integrate() {
    vector<T> result(coefficients.size() + 1, 0);
    for (int i = 0; i < coefficients.size(); ++i) {
        result[i + 1] = coefficients[i] / (i + 1);
    }

    coefficients = result;
    update_coef();
    return *this;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& poly) {
    bool first = true;
    for (int i = poly.power; i >= 0; --i) {
        if (poly.coefficients[i] != 0) {
            if (!first && poly.coefficients[i] > 0) {
                os << "+";
            }
            if (poly.coefficients[i] != 1 || i == 0) {
                os << poly.coefficients[i];
            }
            if (i > 0) {
                os << "x";
                if (i > 1) {
                    os << "^" << i;
                }
            }
            first = false;
        }
    }
    if (first) {
        os << "0";
    }
    return os;
}

template <class T>
istream& operator>>(istream& in, Polynomial<T>& r) {
    cout << "Input power: ";
    in >> r.power;

    r.coefficients.resize(r.power + 1);
    cout << "Input coefficients:" << endl;
    for (int i = 0; i <= r.power; ++i) {
        cout << "Coefficient for x^" << i << ": ";
        in >> r.coefficients[i];
    }
    r.update_coef();
    return in;
}

#endif // !POLYNOM_H