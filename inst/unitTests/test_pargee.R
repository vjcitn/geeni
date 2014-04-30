library(nlme)
library(gee)
library(geepack)
library(HSAUR2)

test_pargee_gaussian <- function()
{
    orthGD <- groupedData(distance ~ age | Subject,
                  data=as.data.frame(Orthodont), 
                  FUN=mean, outer=~ Sex)
    orthGD$Sex <- as.numeric(orthGD$Sex) - 1
    data <- orthGD
#    family <- c(gaussian("identity"), gaussian("log"), gaussian("inverse"))
    formula <- formula(distance ~ age + Sex)

    gee0 <- pargee(orthGD, family=gaussian, binit=c(0, 0, 0))
    gee1 <- gee(formula, data=data, id=Subject, family=gaussian)
    gee2 <- geeglm(formula, data=data, id=Subject, family=gaussian)
    linear <- lm(formula, data=data)

    ## coefficients
    current <- as.vector(round(gee0$coefficients, 4))
    target1 <- unname(round(gee1$coefficients, 4))
    target2 <- unname(round(gee2$coefficients, 4))
    target3 <- unname(round(linear$coefficients, 4))
    checkIdentical(current, target1)
    checkIdentical(current, target2)
    checkIdentical(current, target3)

    ## standard error 
    current <- as.vector(round(sqrt(diag(gee0$robust.variance)), 4))
    target1 <- as.vector(round(sqrt(diag(gee1$robust.variance)), 4)) 
    target2 <- round(summary(gee2)$geese$mean$san.se, 4)
    checkIdentical(current, target1)
    checkIdentical(current, target2)
}

test_pargee_binomial <- function()
{
    ## create new baseline
    data(respiratory)
    resp <- subset(respiratory, month > "0")
    resp$baseline <- 
        rep(subset(respiratory, month == "0")$status, rep(4, 111))
    ## dummy variable for poor respiratory status
    resp$nstat <- as.numeric(resp$status == "good")
    ## remove unused variables for simplicity
    resp <- resp[, !colnames(resp) %in% c("month", "status")] 

    ## centre, treatment, gender and baseline as numeric
    resp$treatment <- as.numeric(factor(resp$treatment)) - 1
    resp$gender <- as.numeric(factor(resp$gender)) - 1
    resp$baseline <- as.numeric(factor(resp$baseline)) - 1
    resp$centre <- as.numeric(resp$centre)
    respGD <- groupedData(nstat ~ centre | subject, 
                          data=resp, FUN=mean, 
                          outer=~ treatment + gender + baseline + age)

    ## coefficients
    formula <- formula(nstat ~ centre + treatment + gender + baseline + age)
    gee0 <- pargee(respGD, family=binomial, binit=rep(0, 6))
    gee1 <- gee(formula, data=resp, family=binomial, id=subject,
                corstr="independence", scale.fix=TRUE, scale.value=1)
    gee2 <- geeglm(formula, data=resp, family=binomial, id=subject)

    current <- as.vector(round(gee0$coefficients, 4))
    current <- current[c(1:4, 6, 5)] ## why is order reversed?
    target1 <- unname(round(gee1$coefficients, 4))
    target2 <- unname(round(gee2$coefficients, 4))
    checkIdentical(current, target1)
    checkIdentical(current, target2)

    ## standard error 
    current <- as.vector(round(sqrt(diag(gee0$robust.variance)), 4))
    current <- current[c(1:4, 6, 5)]
    target1 <- as.vector(round(sqrt(diag(gee1$robust.variance)), 4)) 
    target2 <- round(summary(gee2)$geese$mean$san.se, 4)
    checkIdentical(current, target1)
    checkIdentical(current, target2)
}
