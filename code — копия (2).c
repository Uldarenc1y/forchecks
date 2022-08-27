#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

enum NUMBER_OF_ROOTS
{
    NO_ROOTS, //0
    ONE_ROOT, //1
    TWO_ROOT, //2
    INFINITY_ROOTS = 17, //17 - infinity of roots
};

const double EPS = 1e-8;

struct space // structure of numbers for TestCase
{
    double a, b, c;
    int NROOTS1;
    double x1, x2;
};

int isnul(double i)
{
    return (fabs(i) < EPS); //zero comparison
}

int Parity2Variebles(double a, double b) //check on parity 2 variebles
{
    if (fabs(a - b) < EPS)
        return 1;
    return 0;
}

// docs
//! Solves a square equation ax^2 + bx + c = 0
//! @param [in]    a a-coefficient
//! @param [in]    b b-coefficient
//! @param [in]    c c-coefficient
//! @param [out]  px1 pointer to the 1st root
//! @param [out]  px2 pointer to the 2st root
//! @return NUMBER_OF_ROOTS

int SolutionSquareEq(double a, double b, double c, double* px1, double* px2)
{
    //assert
    assert(isfinite(a) == 1);
    assert(isfinite(b) == 1);
    assert(isfinite(c) == 1);
    assert(!(px1 == 0));
    assert(!(px2 == 0));

    double d = b*b - 4*a*c;

    if (isnul(a)) // a = 0?
        {
        if (isnul(b)) // a = 0, b ?= 0
          {
              if (isnul(c)) // a = 0, b = 0, c ?=0
                return INFINITY_ROOTS;
              return NO_ROOTS;
          }
        if (isnul(c)) // a = 0, c ?= 0
        {
            *px1 = 0;
            return ONE_ROOT;
        }
        else // a = 0
        {
            *px1 = -c/b;
            return ONE_ROOT;
        }

        }
    else if (d < 0) // discriminant ?< 0
        {
         return NO_ROOTS;
         }
    else if (isnul(d)) // discriminant ?= 0
           {
            *px1 = -b/(2*a);
        return ONE_ROOT;
           }
        *px1 = (-b - sqrt(d)) /(2*a); // discriminant ?> 0
        *px2 = (-b + sqrt(d)) /(2*a);

        return TWO_ROOT;

}

void TestCase(struct space Numbers)
{
    double new_x1 = 0, new_x2 = 0; // for comparison with correct number of roots, correct roots
    int  newNROOTS = SolutionSquareEq(Numbers.a, Numbers.b, Numbers.c, &new_x1, &new_x2);
    if (!(newNROOTS == Numbers.NROOTS1))
        printf("Test results are incorrect with a = %lg, b = %lg, c = %lg\nFaked result:Number of roots - %d x1 = %lg x2 = %lg\n\
                  Expected %d number of solutions\nx1 = %lg\nx2 = %lg", Numbers.a, Numbers.b, Numbers.c, newNROOTS, new_x1, new_x2, Numbers.NROOTS1, Numbers.x1, Numbers.x2);
    else if (!(Parity2Variebles(new_x1, Numbers.x1) && Parity2Variebles(new_x2, Numbers.x2)))
        printf("with a = %lg, b = %lg, c = %lg fake x1 = %lg, fake x2 = %lg\n"
                  "expected x1 = %lg, x2 = %lg\n", Numbers.a, Numbers.b, Numbers.c, new_x1, new_x2, Numbers.x1, Numbers.x2);
}

int main()
{
    double a = 0, b = 0, c = 0;
    int p = scanf("%lf%lf%lf", &a, &b, &c);
    while (p < 3)
    {
        while (getchar () != '\n'){};
        printf("enter a, b, c again\n");
        p = scanf("%lf%lf%lf", &a, &b, &c);
    }
    double x1= 0, x2 = 0;
    int n = SolutionSquareEq(a, b, c, &x1, &x2);
    switch(n)
    {
        case NO_ROOTS:
            printf("no roots\n");
            break;

        case ONE_ROOT:
            printf("x = %lg\n", x1);
            break;

        case TWO_ROOT:
            printf("x1 = %lg\nx2 = %lg\n", x1, x2);
            break;

        case INFINITY_ROOTS:
            printf("an infinite number of solutions\n");
            break;
    }
         //for loop on test cases
         //struct

        struct space Numbers [5] = {1, 2, 3, 0,  0, 0,      0, 2, 3, 1, -1.5, 0,     0, 0, 0, 17, 0, 0,      0, 8, 0, 1, 0, 0,   1, 1, 1, 0, 0, 0};
        int number_of_checks = sizeof(Numbers) / sizeof(Numbers[0]);
        for (int i = 0; i < number_of_checks; ++i)
        {
            TestCase(Numbers[i]); // сделать TestCase таким, чтобы он принимал структуры
        }

    return 0;
}

