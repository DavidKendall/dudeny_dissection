## Dudeny Dissection

Julia program to compute and display a Dudeny dissection, following the instructions at
[https://polyhedr.com/dudeny-dissection.html](https://polyhedr.com/dudeny-dissection.html)

A typical execution in the Julia REPL might be:

```julia
julia> include("dudeny.jl)
julia> import .Dudeny as DY
julia> t, r, g, b, y = DY.compute_dissection(5)
julia> DY.show_dissection(t,r,g,b,y)
julia> DY.animate_dissection(t,r,g,b,y)
```
