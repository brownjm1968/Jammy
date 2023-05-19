# Jammy: Julia-Language SAMMY

Simple R-matrix code to demonstrate SAMMY-like cross section calculation in the Julia language.

Notes:
- s-waves only right now
- incident neutrons only right now
- SAMMY is found at: https://code.ornl.gov/RNSD/SAMMY

To install:

```console
git clone https://github.com/brownjm1968/Jammy.git
cd Jammy
julia --project=.
]
instantiate
```

To use: 

```console
import Jammy 
```

If a package dependency is not installed on your system yet (for example `Plots`):

```console
julia
> import Pkg
> Pkg.add("Plots")
```
