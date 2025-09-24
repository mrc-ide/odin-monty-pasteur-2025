# odin-monty-pasteur-2025

Presentation and material for odin-monty introduction to Institut Pasteur 24th September 2025

## Demo script

A demo script containing the code from the slides plus some additional bits is available here: [`odin-monty-demo.R`](odin-monty-demo.R)

### Prerequisites:

* R (4.5.x recommended, 4.3.x or 4.4.x will work)
* RTools on Windows, XCode command line tools on macOS or a functioning C++ toolchain on Linux

Install the packages:

```r
install.packages(
  c("odin2", "dust2", "monty", "decor", "pkgload", "posterior", "brio"),
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org"))
```

Check everything works:

```r
pkgbuild::check_build_tools(debug = TRUE)
```

For more information see the [Installation section of the odin-monty book](https://mrc-ide.github.io/odin-monty/installation.html)

You can now work through the demo script, [`odin-monty-demo.R`](odin-monty-demo.R)
