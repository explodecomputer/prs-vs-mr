# Collider simulation
# Motivated by problem of T2D being adjusted for BMI
# The PRS for adjusted T2D is being tested for assocation with lots of traits
# If T2D -> BMI <- X <- SNP 
# then the genetic score for T2D based on the SNP effects on T2D will associate with X

# To do this simulation properly we need to give T2D some SNPs and 
# evaluate if a PRS threshold in a discovery has an impact


# This hasn't been run yet, but cursory examples show that problem

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


run <- function(nid, nsnp, v1, v2, bxy, bxu, byu)
{
	effs1 <- choose_effects(nsnp, v1)
	effs2 <- choose_effects(nsnp, v2)
	g1 <- make_geno(nid, nsnp, 0.5)
	g2 <- make_geno(nid, nsnp, 0.5)
	x1 <- make_phen(effs1, g1)
	y1 <- make_phen(c(effs2, bxy), cbind(g1, x1))
	u1 <- make_phen(c(bxu, byu), cbind(x1, y1))
	x2 <- make_phen(effs1, g2)
	y2 <- make_phen(c(effs2, bxy), cbind(g2, x2))
	u2 <- make_phen(c(bxu, byu), cbind(x2, y2))

	x1_adj <- residuals(lm(x1 ~ u1))
	x2_adj <- residuals(lm(x2 ~ u2))
	y1_adj <- residuals(lm(y1 ~ u1))
	y2_adj <- residuals(lm(y2 ~ u2))

	# dat1 <- get_effs(x1, y1, g1)
	# dat2 <- get_effs(x2, y2, g2)
	# dat12 <- combine_samples(dat1, dat2) %>% subset(pval.exposure < 5e-8)

	# dat1_adj <- get_effs(x1_adj, y1, g1)
	# dat2_adj <- get_effs(x2_adj, y2, g2)
	# dat12_adj <- combine_samples(dat1_adj, dat2_adj) %>% subset(pval.exposure < 5e-8)

	dat1_adj2 <- get_effs(x1, y1_adj, g1)
	# dat2_adj2 <- get_effs(x2, y2_adj, g2)

	# prs_x2 <- g2 %*% dat1$beta.exposure %>% as.numeric
	# prs_y2 <- g2 %*% dat1$beta.outcome%>% as.numeric
	# prs_x2_adj <- g2 %*% dat1_adj$beta.exposure%>% as.numeric
	# prs_y2_adj <- g2 %*% dat1_adj$beta.outcome%>% as.numeric
	# prs_x2_adj2 <- g2 %*% dat1_adj2$beta.exposure%>% as.numeric
	prs_y2_adj2 <- g2 %*% dat1_adj2$beta.outcome%>% as.numeric

	return(summary(lm(x2 ~ prs_y2_adj2))$coefficients[2,4])

	dat <- data_frame(x1, y1, x2, y2, x1_adj, x2_adj, y1_adj, y2_adj, prs_x2, prs_y2, prs_x2_adj, prs_y2_adj, prs_x2_adj2, prs_y2_adj2)
	return(dat)
}


run(50000, 100, 0.01, 0, 0, sqrt(0.1), sqrt(0.1))


param <- expand.grid(
	nid=50000,
	nsnp=c(10, 100),
	v1 = c(0.01, 0.1),
	v2=0,
	bxy=0,
	bxu=seq(0, 0.1, by=0.01),
	byu=seq(0, 0.1, by=0.01),
	nsim=20
)




with(a, cor(x1, y1))
with(a, cor(x1, y1_adj))
with(a, cor(x1_adj, y1_adj))
with(a, cor(x1_adj, y1))
with(a, cor(prs_x2_adj, y2))
with(a, cor(prs_y2_adj, y2))
with(a, cor(prs_x2_adj, x2))
with(a, cor(prs_x2_adj2, y2))
with(a, cor(prs_y2_adj2, y2))
with(a, cor(prs_x2_adj2, x2))
with(a, cor(prs_y2_adj2, x2))
with(a, cor(prs_y2_adj, y2))

summary(lm(x2 ~ y2_adj, a))

