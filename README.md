# monotone_mixture
### Driver and model definition for a mixture of partially linear models with monotone shape constraints.

To implement: call `flexmix()` from the flexmix package with `mono_reg()` as the model argument. E.g., 
the following call to `flexmix()` produces a model with 6 components, each of which is a partial linear model regressing `Y` on all other 
variables of `df` and with no intercept. The `mon_inc_index` argument to `mono_reg()` instructs the function to estimate a non-parametric, 
monotone-increasing, aka isotonic, function on the second independent variable in `df`.

```R
mod <- flexmix(Y ~ .-1, data = df, k = 6, model = mono_reg(mon_inc_index = 2))
```

flexmix also allows the construction of multiple mixture models within a single object. The following call to `stepFlexmix()` builds 25 mixture models; for each of
`k` equal to 1 through 5 components, `nrep` calls for the construction of 5 models. Now each model contains 2 monotone components -- a monotone-increasing,
or isotonic, relationship between `Y` and the 2nd independent variable of `df`, and a monotone-decreasing, or antitonic, relationship between `Y` and
the 3rd independent variable of `df`.

```R
m2 <- stepFlexmix(Y ~ .-1, data = df, model = mono_reg(mon_inc_index = 2, mon_dec_index = 3), k = 1:5, nrep = 5)
```

For further discussion of the use of the flexmix package, see the guides [here](https://ro.uow.edu.au/cgi/viewcontent.cgi?article=3410&context=commpapers) and
[here](https://cran.rapporter.net/web/packages/flexmix/vignettes/mixture-regressions.pdf), or the following
[blog](https://www.r-bloggers.com/2013/06/estimating-finite-mixture-models-with-flexmix-package/).

### N.B.
The indexing arguments passed to `mono_reg()` are indices of the design matrix constructed by the formula passed to `flexmix`, so they change
based on the design matrix. For example, without an intercept, `mono_inc_index = 2` refers to the 2nd independent variable of the data frame; when 
an intercept *is* included, `mono_inc_index = 2` refers to the 1st independent variable of the data frame.
