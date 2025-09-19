# odin-monty-hpru-2025

Presentation and material for odin-monty introduction to Institut Pasteur 24th September 2025

### Prerequisites:

* R (4.4.x recommended, 4.3.x will work)
* RTools on Windows, XCode command line tools on macOS or a functioning C++ toolchain on Linux

Install the packages:

```r
install.packages(
  c("odin2", "decor", "pkgload", "posterior", "brio"),
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
```

Check everything works:

```r
pkgbuild::check_build_tools(debug = TRUE)
```

For more information see the [Installation section of the odin-monty book](https://mrc-ide.github.io/odin-monty/installation.html)

You can now work through the demo script, [`odin-monty-demo.R`](odin-monty-demo.R)
