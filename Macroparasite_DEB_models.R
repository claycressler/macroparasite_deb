mouseVB <- function(t, y, p) {
    S <- y[1]
    ## max structural length [mm]
    Linf <- unname(p["Linf"])
    ## shape parameter (converts mm to g) [g mm^(-1)]
    alpha <- unname(p["alpha"])
    ## growth rate for host in good condition [day^(-1)]
    gmin <- unname(p["gmin"])

    ## growth rate
    Grate <- 3*gmin*(alpha^(1/3)*Linf*S^(2/3)-S)

    dSdt <- Grate
    return(list(c(dSdt)))
}

## This model is slightly modified to focus on P_2 parasites whose biomass causes a reduction in assimilation efficiency
P2_deb <- function(t, y, p) {
    y <- unname(y)
    ## HOST PARAMETERS
    ## max. SA-specific ingest rate [food g^(-2/3) day^(-1)]
    imax <- unname(p["imax"])
    ## target body condition [dimensionless ratio]
    theta <- unname(p["theta"])
    ## scaling
    eta <- unname(p["eta"])
    # gut transit time [day^(-1)]
    rho <- unname(p["rho"])
    ## biomass conversion efficiencies [dim.]
    epsA <- unname(p["epsA"]) # ingestion to assimilation
    epsG <- unname(p["epsG"]) # assimilate to growth
    epsR <- unname(p["epsR"]) # assimilate to reserve
    epsI <- unname(p["epsI"]) # reserve to immune
    ## max structural length [mm]
    Linf <- unname(p["Linf"])
    ## shape parameter (converts mm to g) [g mm^(-1)]
    alpha <- unname(p["alpha"])
    ## growth rate for host in good condition [day^(-1)]
    gmin <- unname(p["gmin"])
    ## maintenance rates for different tissues [day^(-1)]
    m <- unname(p["m"])
    mc <- unname(p["mc"])
    mi <- unname(p["mi"])
    ## constitutive immune biomass per total weight [dim.]
    k <- unname(p["k"])
    ## induced immune response rate [g^(-1) day^(-1)]
    b <- unname(p["b"])
    ## biomass loss rate of induced immune response [day^(-1)]
    ui <- unname(p["ui"])
    ## host background mortality rate [day^(-1)]
    umin <- unname(p["umin"])
    ## starvation mortality [day^(-1)]
    uacc <- unname(p["uacc"])
    epsAmin <- unname(p["epsAmin"])
    heps <- unname(p["heps"])

    ## PARASITE PARAMETERS
    ## biomass conversion efficiency
    epsP <- unname(p["epsP"])
    ## parasite attack rates [g parasite^(-1) day^(-1)]
    sigmaC <- unname(p["sigmaC"])
    ## parasite handling times [g]
    hC <- unname(p["hC"])
    ## background biomass loss rate [day^(-1)]
    uP <- unname(p["uP"])
    ## biomass loss due to constitutive immune defense [day^(-1) g^(-1)]
    uC <- unname(p["uC"])
    ## biomass loss due to induced immune defense [day^(-1) g^(-1)]
    uI <- unname(p["uI"])
    ## infection age
    tInf <- unname(p["tInf"])

    ## State variables
    G <- y[0]; C <- y[1]; S <- y[2];
    R <- y[3]; Ii <- y[4]; mu <- y[5];
    P2 <- y[6];

    ## assimilation efficiency
    if (t > tInf) {
        epsA <- epsA*(1-epsAmin*P2/(heps+P2))
    } else epsA <- epsA
    ## ingestion rate
    In <- imax*S^(2/3)/(1+exp(eta*(R/S-theta)))
    ## growth rate
    Grate <- 3*gmin*(alpha^(1/3)*Linf*S^(2/3)-S)
    ## cost of growth
    CG <- epsG*Grate
    ## constitutive immune defense
    Ic <- k*(S+R)
    ## maintenance rate
    M <- m*(S+R) + mc*Ic + mi*Ii

    ## Rate equations
    dGdt <- In - rho*G
    dCdt <- rho*(1-epsA)*G - rho*C
    dSdt <- Grate
    dRdt <- epsR*(rho*epsA*G-M-CG)
    dmudt <- umin+uacc*max(theta*S/R-1,0)
    dIidt <- 0
    dP2dt <- 0
    if (t > tInf) {
        dCdt <- dCdt - sigmaC*P2*C/(hC+C)
        dRdt <- dRdt - b*R*P2
        dIidt <- dIidt + epsI*b*R*P2 - ui*Ii
        dP2dt <- dP2dt + epsP*sigmaC*P2*C/(hC+C) - uP*P2 - uC*Ic*P2 - uI*Ii*P2
    }
    return(list(c(G=dGdt, C=dCdt, S=dSdt, R=dRdt, Ii=dIidt, mu=dmudt, P2=dP2dt), ingest=In, assim=epsA, immalloc=b*R*P2))

}
