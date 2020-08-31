#' Probability sub-sampling procedure
#'
#' This procedure uses reference information at the target population level to
#' select a sub-sample which would be similar to a random selection from it
#'
#' @param Input.Sample A data frame with the initial sample that needs subs-
#'   sampling. It needs to have the columns: Code, Region, Beds, QS.
#' @param Reference.Data Target population level reference data, with
#'   information on number of Beds and Region of all hospitals.
#' @param n.required  The size of the final sub-sample.
#' @param n.quantiles Number of quantiles into which categorize continuous
#'   variables, like the hospital size.
#' @param QS.weight Scaling factor of the influence of the QS in the algorithm.
#'   Set to 0 to remove the effect of the QS.
#' @param method Whether to select the hospitals in a deterministic way after
#'   arranging them by their score, or in a probabilistic way using the score as
#'   a sampling weight.
#'
#' @import dplyr
#' @import Hmisc
#' @import magrittr
#' @import scales
#'
#' @return A subset of the data frame passed in Input.Sample with n.required
#'   rows
#' @export
#'
#' @examples
#'
#' # Create a subsample of 55 hospitals
#' subsample.probability(Sample.Data, Reference.Data, n.required = 55)
#'
#' # Set QS.weight to 0 to not use the QS in the sampling
#' subsample.uniform(Sample.Data, Reference.Data, n.required = 55, QS.weight = 0)
#'
#' # Create a subsample of 55 hospitals with weighted random selection
#' subsample.probability(
#'     Sample.Data, Reference.Data, n.required = 55, method = 'random'
#' )

probability.sampling <- function(Input.Sample, Reference.Data, n.required,
																	n.quantiles = 10, QS.weight = 1,
																	method = c('arrange', 'random')){

	library(dplyr)
	library(Hmisc)
	library(magrittr)
	library(scales)

	method <- match.arg(method)

	# Definition of quantiles in the distribution of number of beds according to
	# reference data
	quantiles <- quantile(Reference.Data$Beds, seq(0, 1, length.out = n.quantiles)) %>%
		round()

	# Definition of the distributional blocks at the target population level
	P_country <- Reference.Data %>%
		count(Block = Hmisc::cut2(Beds, quantiles) %>% paste('-', Region)) %>%
		mutate(Prob = n / sum(n)) %>%
		with(magrittr::set_names(Prob, Block))

	Selection <- Input.Sample %>%
		mutate(
			# Identification of the country level blocks in the sample
			Block = Hmisc::cut2(Beds, quantiles) %>% paste('-', Region),
			# Association of P_country to the hospital in the sample
			Prob = P_country[Block],
			# Creation of the quality weight after rescaling of the QS
			QS.rescale = 1 - scales::rescale(QS),
			# Definition of the final score
			Score = Prob * QS.rescale^QS.weight
		) %>%
		filter(!is.na(Score)) # Remove blocks that do not appear in the national list

	if (method == 'arrange') {
		Selection %>%
			arrange(desc(Score)) %>%
			head(n.required)
	} else {
		slice_sample(Selection, n = n.required, weight_by = Score)
	}
}