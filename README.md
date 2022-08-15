
## Dynamic risk prediction triggered by intermediate events using survival tree ensembles

This repository contains example codes used to run the simulation in the
manuscript “Dynamic risk prediction triggered by intermediate events
using survival tree ensembles” described in
<https://arxiv.org/abs/2011.07175>.

This version was built in March, 2022. The contents of the repository
are as follows.

-   The `/codes` folder contains all the required `R` scripts to
    reproduce results in our examples.
-   The `/Examples` folder contains `Rmd` documents that we used to
    built online vignettes whose links are attached below.

## Example

Depending on the landmark time and how the longitudinal predictors are
measured, we provide the following vignettes to introduce the
implementations of the proposed ensemble methods.

Simulation settings and essential functions:

-   [Functions to set up
    simulation](https://htmlpreview.github.io/?https://github.com/stc04003/DynmaicRisk/blob/main/Examples/data.html).
-   [Functions for ensemble
    procedures](https://htmlpreview.github.io/?https://github.com/stc04003/DynmaicRisk/blob/main/Examples/addon.html).

Examples:

-   [Random landmark time with longitudinal predictors measured at 1, 2,
    …, 5.](https://htmlpreview.github.io/?https://github.com/stc04003/DynmaicRisk/blob/main/Examples/sim3A.html).
-   [Fixed landmark time with longitudinal predictors measured at the
    intermediate
    event.](https://htmlpreview.github.io/?https://github.com/stc04003/DynmaicRisk/blob/main/Examples/sim3B.html).
-   [Random landmark time with longitudinal predictors measured at the
    intermediate
    event.](https://htmlpreview.github.io/?https://github.com/stc04003/DynmaicRisk/blob/main/Examples/sim3C.html).
