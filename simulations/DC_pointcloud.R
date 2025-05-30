# -------------------------------------------------------------------------
# Simulation Study: Point cloud data with intrinsic distance
#
# Manuscript reference : Section 5.2
# Figure reproduced    : Figure 5
# Table reproduced     : Table 1
# -------------------------------------------------------------------------


## Import main functions
source("src/DC_mainfunctions.R")
source("simulations/DC_simulations_datagen.R")


#########################################
## Point Cloud Data B1: Positively curved manifold
k = 1000 # Number of samples 

## Generate simulation data from point cloud on positively curved manifold
df.pc.pos = dat.gen.pc.pos(k)

## Estimate intrinsic curvature and test statistics (main function)
curv.fit.pc.pos = intrinsic.curv.est(df.pc.pos, L = 4)
cov.est.pc.pos = curv.fit.pc.pos$cov.normalized
Vm.est.pc.pos = curv.fit.pc.pos$Vm
Vf.est.pc.pos = curv.fit.pc.pos$Vf

## Confidence intervals for the intrinsic curvatures rho_{I}
# alpha = 0.01 (0.072, 0.132)
C.I.lower.pc.pos.01 = curv.fit.pc.pos$rho - qnorm(0.995)*curv.fit.pc.pos$sd/sqrt(k)
C.I.upper.pc.pos.01 = curv.fit.pc.pos$rho + qnorm(0.995)*curv.fit.pc.pos$sd/sqrt(k)

# alpha = 0.05 (0.080, 0.125)
C.I.lower.pc.pos.05 = curv.fit.pc.pos$rho - qnorm(0.975)*curv.fit.pc.pos$sd/sqrt(k)
C.I.upper.pc.pos.05 = curv.fit.pc.pos$rho + qnorm(0.975)*curv.fit.pc.pos$sd/sqrt(k)

# alpha = 0.1 (0.083, 0.122)
C.I.lower.pc.pos.1 = curv.fit.pc.pos$rho - qnorm(0.95)*curv.fit.pc.pos$sd/sqrt(k)
C.I.upper.pc.pos.1 = curv.fit.pc.pos$rho + qnorm(0.95)*curv.fit.pc.pos$sd/sqrt(k)

## Table 1
ci.tbl.pos <- data.frame(
  Level = c("99%", "95%", "90%"),
  Lower = c(C.I.lower.pc.pos.01,
            C.I.lower.pc.pos.05,
            C.I.lower.pc.pos.1),
  Upper = c(C.I.upper.pc.pos.01,
            C.I.upper.pc.pos.05,
            C.I.upper.pc.pos.1)
)

ci.tbl.pos <- transform(ci.tbl.pos,
                    Lower = round(Lower, 3),
                    Upper = round(Upper, 3))

print("Point Cloud Data B1: Positively curved manifold")
print(ci.tbl.pos, row.names = FALSE)


#########################################
## Point Cloud Data B2: Negatively curved manifold

## Generate simulation data from point cloud on positively curved manifold
df.pc.neg = dat.gen.pc.neg(k)

## Estimate intrinsic curvature and test statistics (main function)
curv.fit.pc.neg = intrinsic.curv.est(df.pc.neg, L = 4)
cov.est.pc.neg = curv.fit.pc.neg$cov.normalized
Vm.est.pc.neg = curv.fit.pc.neg$Vm
Vf.est.pc.neg = curv.fit.pc.neg$Vf

## Confidence intervals for the intrinsic curvatures rho_{I}
# alpha = 0.01 (-0.098, -0.072)
C.I.lower.pc.neg.01 = curv.fit.pc.neg$rho - qnorm(0.995)*curv.fit.pc.neg$sd/sqrt(k)
C.I.upper.pc.neg.01 = curv.fit.pc.neg$rho + qnorm(0.995)*curv.fit.pc.neg$sd/sqrt(k)

# alpha = 0.05 (-0.095, -0.075)
C.I.lower.pc.neg.05 = curv.fit.pc.neg$rho - qnorm(0.975)*curv.fit.pc.neg$sd/sqrt(k)
C.I.upper.pc.neg.05 = curv.fit.pc.neg$rho + qnorm(0.975)*curv.fit.pc.neg$sd/sqrt(k)

# alpha = 0.1 (-0.093, -0.077)
C.I.lower.pc.neg.1 = curv.fit.pc.neg$rho - qnorm(0.95)*curv.fit.pc.neg$sd/sqrt(k)
C.I.upper.pc.neg.1 = curv.fit.pc.neg$rho + qnorm(0.95)*curv.fit.pc.neg$sd/sqrt(k)

## Table 1
ci.tbl.neg <- data.frame(
  Level = c("99%", "95%", "90%"),
  Lower = c(C.I.lower.pc.neg.01,
            C.I.lower.pc.neg.05,
            C.I.lower.pc.neg.1),
  Upper = c(C.I.upper.pc.neg.01,
            C.I.upper.pc.neg.05,
            C.I.upper.pc.neg.1)
)

ci.tbl.neg <- transform(ci.tbl.neg,
                    Lower = round(Lower, 3),
                    Upper = round(Upper, 3))

print("Point Cloud Data B2: Negatively curved manifold")
print(ci.tbl.neg, row.names = FALSE)

#########################################
## Point Cloud Data B3: Flat manifold

## Generate simulation data from point cloud on positively curved manifold
df.pc.flat = dat.gen.pc.flat(k)

## Estimate intrinsic curvature and test statistics (main function)
curv.fit.pc.flat = intrinsic.curv.est(df.pc.flat, L = 4)

## Plot confidence region for the joint distribution of intrinsic metric and Frechet variance
cov.est.pc.flat = curv.fit.pc.flat$cov.normalized
Vm.est.pc.flat = curv.fit.pc.flat$Vm
Vf.est.pc.flat = curv.fit.pc.flat$Vf

## Confidence intervals for the intrinsic curvatures rho_{I}
# alpha = 0.01 (-0.023, 0.051)
C.I.lower.pc.flat.01 = curv.fit.pc.flat$rho - qnorm(0.995)*curv.fit.pc.flat$sd/sqrt(k)
C.I.upper.pc.flat.01 = curv.fit.pc.flat$rho + qnorm(0.995)*curv.fit.pc.flat$sd/sqrt(k)

# alpha = 0.05 (-0.014, 0.042)
C.I.lower.pc.flat.05 = curv.fit.pc.flat$rho - qnorm(0.975)*curv.fit.pc.flat$sd/sqrt(k)
C.I.upper.pc.flat.05 = curv.fit.pc.flat$rho + qnorm(0.975)*curv.fit.pc.flat$sd/sqrt(k)

# alpha = 0.1 (-0.009, 0.038)
C.I.lower.pc.flat.1 = curv.fit.pc.flat$rho - qnorm(0.95)*curv.fit.pc.flat$sd/sqrt(k)
C.I.upper.pc.flat.1 = curv.fit.pc.flat$rho + qnorm(0.95)*curv.fit.pc.flat$sd/sqrt(k)

## Table 1
ci.tbl.flat <- data.frame(
  Level = c("99%", "95%", "90%"),
  Lower = c(C.I.lower.pc.flat.01,
            C.I.lower.pc.flat.05,
            C.I.lower.pc.flat.1),
  Upper = c(C.I.upper.pc.flat.01,
            C.I.upper.pc.flat.05,
            C.I.upper.pc.flat.1)
)

ci.tbl.flat <- transform(ci.tbl.flat,
                        Lower = round(Lower, 3),
                        Upper = round(Upper, 3))

print("Point Cloud Data B3: Flat manifold")
print(ci.tbl.flat, row.names = FALSE)

