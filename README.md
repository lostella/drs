## Douglas-Rachford splitting

This repository contains implementations of the Douglas-Rachford operators
splitting method for solving convex optimization problems of the form

```
minimize f(x) + g(x)
```

The implementations and the proposed tests reflect the analysis contained in
[1]. In particular, a *fast* version of the method is provided for problems
where f(x) is quadratic.

### References

[1] P. Patrinos, L. Stella, A. Bemporad, “Douglas-Rachford Splitting: Complexity
Estimates and Accelerated Variants,” [arXiv:1407.6723](http://arxiv.org/abs/1407.6723) [math.OC], Sep. 2014.
