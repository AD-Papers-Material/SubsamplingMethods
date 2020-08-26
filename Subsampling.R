

#' Uniform sub-sampling procedure
#'
#' Create a sub-sample which tends to be have a uniform distribution with
#' respect to hospital location and hospital size as number of beds.
#'
#' @param Input.Sample A data frame with the initial sample that needs subs-
#' sampling. It needs to have the columns: Code, Region, Beds, QS (if required).
#' @param n.required  The size of the final sub-sample.
#' @param n.quantiles Number of quantiles into which categorize continuous
#' variables, like the hospital size.
#' @param use.QS Whether to choose the sequence of hospital to evaluate for
#' selection based on the quality of their data or randomly.
#'
#' @import dplyr
#' @import Hmisc
#'
#' @return A subset of the data frame passed in Input.Sample with n.required rows.
#' @export
#'
#' @examples
#'
#' # Create a subsample of 56 hospitals using the QS
#' subsample.uniform(Sample.Data, n.required = 56)
#'
#' # The same but this time hospitals are chosen randomly
#' subsample.uniform(Sample.Data, n.required = 56)

subsample.uniform <- function(Input.Sample, n.required, n.quantiles = 4, use.QS = T){
	library(dplyr)
	library(Hmisc)

	Hospitals <- Input.Sample %>%
		transmute(Code, Region, Beds = Hmisc::cut2(Beds, g = n.quantiles), QS = if (use.QS) QS else 1) # Prepare data by discretizing continuous variables like the number of acute beds and by changing QS to a fixed value if not to be used

	Selected.hospitals <- c()
	Candidates <- data.frame()

	for (i in 1:n.required) { # Until n.required is reached..
		if (nrow(Candidates) == 0) { # If the candidate list is empty, rebuilt it from the non-selected hospitals
			Candidates <- Hospitals %>%
				filter(!(Code %in% Selected.hospitals)) %>% # Remove already selected hospitals
				sample_frac() %>% # Permute order, useful only if QS is not used
				arrange(QS)
		}

		# Extract the first hospital of the temporary list and add it to the list of selected hospitals
		Extracted.hospital <- Candidates[1,]
		Selected.hospitals <- c(Selected.hospitals, Extracted.hospital$Code)

		# Remove from the temporary list all hospitals in the same location/size block of the extracted hospital
		Candidates <- Candidates %>%
			filter(!(Region %in% Extracted.hospital$Region & Beds %in% Extracted.hospital$Beds))
	}

	Input.Sample %>% filter(Code %in% Selected.hospitals) # Filter the initial data by the selected hospital codes
}



#' Probability sub-sampling procedure
#'
#' This procedure uses reference information at the target population level to
#' select a sub-sample which would be similar to a random selection from it
#'
#' @param Input.Sample A data frame with the initial sample that needs subs-
#' sampling. It needs to have the columns: Code, Region, Beds, QS.
#' @param Reference.Data Target population level reference data, with
#' information on number of Beds and Region of all hospitals.
#' @param n.required  The size of the final sub-sample.
#' @param n.quantiles Number of quantiles into which categorize continuous
#' variables, like the hospital size.
#' @param QS.weight Scaling factor of the influence of the QS in the algorithm.
#' Set to 0 to remove the effect of the QS.
#' @param method Whether to select the hospitals in a deterministic way after
#' arranging them by their score, or probabilistically using the score as a
#' sampling weight.
#'
#' @import dplyr
#' @import Hmisc
#' @import magrittr
#' @import scales
#'
#' @return A subset of the data frame passed in Input.Sample with n.required rows
#' @export
#'
#' @examples
#'
#' # Create a subsample of 56 hospitals
#' subsample.probability(Sample.Data, Reference.Data, n.required = 56)
#'
#' # Set QS.weight to 0 to not use the QS in the sampling
#' subsample.uniform(Sample.Data, Reference.Data, n.required = 56, QS.weight = 0)
#'
#' # Create a subsample of 56 hospitals with weighted random selection
#' subsample.probability(Sample.Data, Reference.Data, n.required = 56, method = 'random')

subsample.probability <- function(Input.Sample, Reference.Data, n.required, n.quantiles = 10, QS.weight = 1, method = c('arrange', 'random')){

	library(dplyr)
	library(Hmisc)
	library(magrittr)
	library(scales)

	method <- match.arg(method)

	# Definition of quantiles in the distribution of number of beds according to reference data
	quantiles <- quantile(Reference.Data$Beds, seq(0, 1, length.out = n.quantiles)) %>% round

	P_country <- Reference.Data %>%
		count(Block = Hmisc::cut2(Beds, quantiles) %>% paste('-', Region)) %>%
		mutate(Prob = n / sum(n)) %>%
		with(magrittr::set_names(Prob, Block))

	Selection <- Input.Sample %>%
		mutate(
			Block = Hmisc::cut2(Beds, quantiles) %>% paste('-', Region), # Identification of the country level blocks in the sample
			Prob = P_country[Block], # Assoction of P_country to the hospital in the sample
			QS.rescale = 1 - scales::rescale(QS), # Creation of the quality weight after rescaling of the QS
			Score = Prob * QS.rescale^QS.weight # Definition of the final score
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