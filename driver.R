# imports

#TODO general formula is passed to all mods, formula = form passed to individual model construction ?
m1 <- flexmix(y ~ 1, k = 4, model = mono_reg(formula = form), data = dat)

# stepFlexmix runs flexmix for nrep=n iterations.
m7 <- stepFlexmix(yp ~ x + I(x^2), data = NPreg,
                  +    control = list(iter.max=1,
                                      minprior=0.1,
                                      tolerance=0.001,
                                      verbose=0,
                                      classify="weighted",
                                      nrep=5), k = 1:5, nrep = 5)