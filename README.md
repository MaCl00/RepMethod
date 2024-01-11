# Documentation - RepMethod
 An Automated Repertoire Method for solving linear recurrence relations
 Requires MatLab and the Symbolic Math Toolbox Add-On.
## Syntax
result = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree, optional_flags)
[C, N] = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree)
[__] = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree, true)
[__] = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree, true, false)
[__] = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree, false, false)
## Description
result = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree) returns the solutions to the recurrence relation implied by substitute and nonHomogeneous which are present in the span of functions
[C, N] = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree) saves the coefficients of the solutions in the matrix C and the names of the functions in N
[__] = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree, true) gives verbose information during the execution. The default is false.
[__] = repertoire(functions, @substitute, nonHomogeneous, precision, start, rec_degree, poly_degree, false, false) does not quit searching for solutions if the probability for further solutions is low. The default value is true (the search **is** quit when the probability for further solutions is low)
## Examples
For detailed examples, refer to the directory *examples_commented*
## Input Arguments
### functions - function guesses
A list of symbolic functions that the user suspects as solutions. Products are also checked (e.g. if *f* and *g* are in *functions*, *f* * *g* will be added to *functions*). The program can find solutions within the span of *functions*. Be sure to use the symbolic variable *x*.
### @substitute - homogeneous part of the linear recurrence relation
The function implies the homogeneous part of the linear recurrence relation. E.g. if $a_n = n\cdot a_{n-1} + 2^n \cdot a_{n-2} + n!$ is the recurrence, then substitute should be 

```
function ret = substitute(n, func_val)
   ret = func_val(:,1) - sym(n) * func_val(:,2) + 2^sym(n) * func_val(:,3);
end
```

### nonHomogeneous - symbolic function of the non-homogeneous part of the recurrence relation
The non-homogeneous part of the recurrence relation. E.g. if $a_n = n\cdot a_{n-1} + 2^n \cdot a_{n-2} + n!$ is the recurrence, then *nonHomogeneous* is *factorial(sym(x))* (be sure to define beforehand *syms x*). The non-homogeneous solution will be saved in the first entry of *C*, if no non-homogeneous solution is found, the first entry is *0*. In case of a homogeneous recurrence, set *nonHomogeneous* to $0$.
### precision - calculation precision
*precision* specifies the precision of all calculations (in significant decimal figures). Usually, a precision of ~400 yields good results. 
### start - first integer at which functions get evaluated
*start* is the first integer at which all functions in *functions* are evaluated. Subsequent points are the consecutive integers starting from *start*.
### rec_degree - denotes the *order*$+1$ of the recurrence relation
In the example $a_n = n\cdot a_{n-1} + 2^n \cdot a_{n-2} + n!$ *rec_degree* would be $3$. It is equivalent to the number of terms in the second parameter of substitute. If a recurrence relation references (nearly) all previous terms, for example
$$a_n=\sum_{i=0}^{n-1}  a_i$$
then *rec_degree* should be $-1$.
### poly_degree - denotes the number of the polynomials that are guessed
*poly_degree* specifies how many polynomials are guessed. For example, if *poly_degree* is $3$ the program adds $1,x$ and $x^2$ to *functions*. Polynomials are also pairwise multiplied with other functions in *functions*. For example, if the user provides *functions* with two functions $f(x),g(x)$ and sets *poly_degree* to $3$, the program guesses: $f(x),g(x),f(x)\cdot g(x), x\cdot f(x),x\cdot g(x), x\cdot f(x)\cdot g(x),x^2\cdot f(x),x^2\cdot g(x),x^2\cdot f(x)\cdot g(x),1,x,x^2$
### optionalFlag1=false - verbose flag
An optional flag that specifies whether the program provides information in the console. By default, the program does not provide additional information in the console.
### optionalFlag2=true - search quit flag
An optional flag that specifies whether the search is quit, if the program calculates a low probability for further solutions. By default the search is quit, if the probability for further solutions is low.
