# ExomeDepthExample

Example for running ExomeDepth WES samples, accounting for relatedness in the reference sample distribution

1. Clone this repository by running `git clone https://github.com/PankratzLab/ExomeDepthExample.git`
1. Install ExomeDepth
	- https://github.com/vplagnol/ExomeDepth
	- `install.packages("ExomeDepth")`
	- read the docs https://cran.r-project.org/web/packages/ExomeDepth/index.html .... but it looks like the package was removed from cran
```
Package ‘ExomeDepth’ was removed from the CRAN repository.

Removed on 2022-10-09 for misrepresentation of authorship and copyright holders.

Please use the canonical form https://CRAN.R-project.org/package=ExomeDepth to link to this page.

```

3. Install openxlsx
	- `install.packages("openxlsx")`
4. Adopt exomeDepthRelated.R to correct paths, and then run the script
