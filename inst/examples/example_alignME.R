\dontrun{

data(MAPK)

## Run alignME with standard arguments
out <- alignME(data = MAPK,
               model = "ys/sj",
               errmodel = "sigmaR*value",
               fixed = ys~Condition,
               latent = sj~Experiment,
               error = sigmaR~1,
               log = TRUE)
plot1(out)
plot2(out)
plot3(out)

## Assume equal variance on all gels
out <- alignME(data = MAPK,
               model = "ys/sj",
               errmodel = "sigma0",
               fixed = ys~Condition,
               latent = sj~Experiment,
               error = sigma0~1,
               log = TRUE)
plot2(out)

## Estimate with offset
out <- alignME(data = MAPK,
               model = "ys/sj+bj",
               errmodel = "sigmaR*value",
               fixed = ys~Condition,
               latent = sj+bj~Experiment,
               error = sigmaR~1,
               log = TRUE)
plot2(out)

## Align data on log-scale
logMAPK <- MAPK
logMAPK$value <- log(MAPK$value)
out <- alignME(data = logMAPK,
               model = "log(ys/sj)",
               errmodel = "sigmaR",
               fixed = ys~Condition,
               latent = sj~Experiment,
               error = sigmaR~1,
               log = TRUE)
plot1(out)

## Align data on log-scale with mixed error model
logMAPK <- MAPK
logMAPK$value <- log(MAPK$value)
out <- alignME(data = logMAPK,
               model = "log(ys/sj)",
               errmodel = "sigmaR + sigma0*exp(-value)",
               fixed = ys~Condition,
               latent = sj~Experiment,
               error = sigmaR+sigma0~1,
               log = TRUE)
plot1(out)

}
