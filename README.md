## Dudeny Dissection

Julia program to compute and display a Dudeny dissection, following the instructions at
[https://polyhedr.com/dudeny-dissection.html](https://polyhedr.com/dudeny-dissection.html)

### Installation
Clone the repository, enter it and activate the project, e.g.

```shell
$ git clone https://github.com/DavidKendall/dudeny_dissection
$ cd dudeny-dissection
$ julia --project=.
```

### Usage
A typical execution in the Julia REPL might be:

```julia-repl
julia> include("dudeny.jl)
julia> import .Dudeny as DY
julia> t, r, g, b, y = DY.compute_dissection(5)
julia> DY.show_dissection(t,r,g,b,y)
julia> DY.animate_dissection(t,r,g,b,y)
```

### Notes
1. The first time you run the project there will be some delay while Julia installs
   the packages that are required.
