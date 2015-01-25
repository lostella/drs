## Douglas-Rachford splitting

This repository contains implementations of the Douglas-Rachford operator
splitting method for solving convex optimization problems of the form

```
minimize f(x) + g(x)
```

The proposed implementations and tests reflect the analysis contained in
[1]. In particular, a *fast* version of the method is provided for problems
where f(x) is quadratic.

The problems for which the algorithms are implemented are:
* Box constrained QPs
* L1-regularized least square regression

### References

[1] P. Patrinos, L. Stella, A. Bemporad, “Douglas-Rachford Splitting: Complexity
Estimates and Accelerated Variants,” [arXiv:1407.6723](http://arxiv.org/abs/1407.6723) [math.OC], Sep. 2014.
In *Proceedings of the 53nd IEEE Conference on Decision and Control*.
