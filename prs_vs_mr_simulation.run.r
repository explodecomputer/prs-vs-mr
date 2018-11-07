library(simulateGP)
library(TwoSampleMR)

combine_samples <- function(dat1, dat2)
{
	a <- subset(dat1, select=c("SNP", grep("exposure", names(dat1), value=TRUE)))
	b <- subset(dat2, select=c("mr_keep", grep("outcome", names(dat1), value=TRUE)))
	return(cbind(a,b))
}

prs <- function(dat1, g2, y2)
{
	score <- g2 %*% dat1$beta.exposure
	sc <- summary(lm(y2 ~ score))
	return(data.frame(
		beta = coefficients(sc)[2,1],
		se = coefficients(sc)[2,2],
		pval = coefficients(sc)[2,4]
	))
}


run <- function(nid, nsnp, v1, v2, bxy)
{
	effs1 <- choose_effects(nsnp, v1)
	effs2 <- choose_effects(nsnp, v2)
	g1 <- make_geno(nid, nsnp, 0.5)
	g2 <- make_geno(nid, nsnp, 0.5)
	x1 <- make_phen(effs1, g1)
	y1 <- make_phen(c(effs2, bxy), cbind(g1, x1))
	x2 <- make_phen(effs1, g2)
	y2 <- make_phen(c(effs2, bxy), cbind(g2, x2))

	dat1 <- get_effs(x1, y1, g1)
	dat2 <- get_effs(x2, y2, g2)
	dat12 <- combine_samples(dat1, dat2) %>% subset(pval.exposure < 5e-8)
	if(nrow(dat12) == 0)
	{
		a <- list()
		a$pval <- NA
	} else {
		a <- mr(dat12, method_list="mr_ivw")
	}
	ind <- dat1$pval.exposure < 1e-5
	if(sum(ind) == 0)
	{
		b <- list()
		b$pval <- NA
	} else {
		b <- prs(dat1[ind,], g2[,ind], y2)
	}
	return(list(mr=a$pval, prs=b$pval))
}

param <- expand.grid(
	v1=c(0.15),
	# v2=c(0, seq(0.01, 0.2, by=0.01)),
	bxy=c(0, 0.4),
	nsim=1:1000,
	v2 = NA,
	mr = NA,
	prs = NA
)
param$v2 <- runif(nrow(param), min=0, max=0.2)
for(i in 1:nrow(param))
{
	message(i, " of ", nrow(param))
	out <- run(10000, 50, param$v1[i], param$v2[i], param$bxy[i])
	param$mr[i] <- out$mr
	param$prs[i] <- out$prs
}
save(param, file="prs_vs_mr_simulation.rdata")
