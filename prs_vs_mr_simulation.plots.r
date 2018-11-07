library(tidyverse)

load("prs_vs_mr_simulation.rdata")
paraml <- gather(param, "key", "value", mr, prs)
paraml$lab <- "Null model"
paraml$lab[paraml$bxy == 0.4] <- "Causal model"
out <- paraml %>% mutate(plei_bin=cut(v2/v1, 50)) %>%
	group_by(plei_bin, lab, bxy, key) %>%
	summarise(n=n(), prop=sum(value < 1e-10)/n)
paraml$value[paraml$lab == "Null model" & paraml$value < 1e-20] <- 1e-20

ggplot(paraml, aes(x=v2 / v1, y=-log10(value))) +
geom_point(aes(colour=key), size=0.2) +
geom_smooth(aes(colour=key)) +
facet_grid(lab ~ ., scale="free_y") +
labs(x="Proportion of genetic variation due to balanced horizontal pleiotropy", y="-log10 p-value", colour="method")
ggsave("prs_vs_mr_simulation.plot1.pdf")


ggplot(out, aes(x=plei_bin, y=prop)) +
geom_point(aes(colour=key)) +
geom_smooth(aes(colour=key)) +
facet_grid(lab ~ ., scale="free_y") +
labs(x="Proportion of genetic variation due to balanced horizontal pleiotropy", y="-log10 p-value", colour="method")
ggsave("prs_vs_mr_simulation.plot2.pdf")

