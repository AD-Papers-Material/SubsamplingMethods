#' Title
#'
#' @param Sample.Data
#' @param Reference.Data
#' @param n_groups
#' @param to.select
#' @param qs.weight
#'
#' @return
#' @export
#'
#' @examples
probability.sampling <- function(Sample.Data, Reference.Data, n_groups = 10, to.select, qs.weight = 1) {

	# Definition of quantiles in the distribution of number of beds according to the n_groups parameter
	quantiles <- quantile(Refence.Data$Beds, seq(0, 1, length.out = n_groups)) %>% round

	#
	Prob.groups <- Reference.Data %>%
		with(
			cut2(Beds, quantiles) %>% paste('-', Region)
		) %>% {table(Group = .)} %>% as.data.frame() %>% mutate(Prob = Freq / sum(Freq))

	Prob.groups <- Prob.groups$Prob.country %>% set_names(Prob.groups$Group)

	Sample.Data %>%
		filter(!is.na(QS)) %>%
		mutate(
			P.group = cut2(Beds, quantiles) %>% paste('-', Region),
			Prob = Prob.groups[P.group],
			QS.rescale = 1 - rescale(QS),
			Prob.score = Prob * QS.rescale^qs.weight
		) %>% sample_frac() %>%
		arrange(desc(Prob.score)) %>%
		head(to.select)
}