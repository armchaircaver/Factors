C / C++ program suite to factorise 64 bit numbers as quickly as possible.

The program suite incorporates a number of techniques to produce speedy factorisation of numbers, including:

- Using libdivide to speed up the trial division part of the factorisation.
    Libdivide structures for small primes are pre-calculated in the initialisation, and used to test 
    divisibility using "division by multiplication". See https://libdivide.com/ for details.
    
- Implementing modular multipication speedily for moduli below 2^63 using the Rutten-Eekelen technique  
https://www.researchgate.net/publication/220493198_Efficient_and_Formally_Proven_Reduction_of_Large_Integers_by_Small_Moduli
Rutten-Eekelen modular multiplication is used in the Miller Rabin tests.

- Using the Fast Primality Testing technique developed by Michal Foriˇsek and Jakub Janˇcina
http://ceur-ws.org/Vol-1326/020-Forisek.pdf to perform Miller-Rabin primality testing.
Their 64 bit implementation uses Montgomery Multiplication, but my implementation of Miller-Rabin only uses this for larger numbers above 2^43, and uses Rutten-Eekelen
for smaller numbers.

- Using Montgomery Multiplication in the Pollard Rho part of the factorisatsion.



    
    

