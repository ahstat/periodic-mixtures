# Mélanges périodiques

README


- f_g_Ff_Fg.R: direct computation of f, g, Ff and Fg
- the related tests give what is into the plot1 folder (Base functions section)

- Sf_Sg.R: sum of f and g computed with direct approx, Fourier approx, and closed-form
- the related tests give what is into the plot2 folder (Sum of the base functions section), and also
  the step-by-step for the closed form (Closed-form expressions section)

- elements.R: 
  build_f_g_Ff_Fg(), to get a list containing the functions f,g,Ff and Fg for a certain type and sigma
  and build_Sf_Sg to get in addition the computations of the functions Sf and Sg

- helpers.R: some function, and also some numerical root computations

- helpers_for_tests.R: for automazation of loop of tests, along with save plots to png

- dynamics.R: ok until gaussian case

TODO: numerical results for gaussian and derivatives + plot + explained
step0: verifier code et faire tous les plots
step05: verifier que les 3 normalisations choisies sont correctes
step1: faire numerical
step2: faire closed form if easy
  - gaussian normalization 1
  - gaussian normalization 2
  - gaussian normalization 3



### Quick compilation command

``` r
rm(list = ls())
devtools::document()
devtools::test()
devtools::build()
library(periodicmixtures)

# Compile the documentation  'Ctrl + Shift + D'
# Build and Reload Package:  'Ctrl + Shift + B'
# Check Package:             'Ctrl + Shift + E'
# Test Package:              'Ctrl + Shift + T'

https://ahstat.github.io/Periodic-mixtures/?version=1234567
```
