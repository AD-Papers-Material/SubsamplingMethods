

#' Uniform sub-sampling procedure
#'
#' Create a sub-sample which tends to be have a uniform distribution with
#' respect to hospital location and hospital size as number of beds.
#'
#' @param InputSample A data frame with the initial sample that needs subs-
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
#' @return A subset of the data frame passed in InputSample.
#' @export
#'
#' @examples
#'
#' # Create a subsample of 56 hospitals using the QS
#' subsample.uniform(SampleData, n.required = 56)
#'
#' # The same but this time hospitals are chosen randomly
#' subsample.uniform(SampleData, n.required = 56)

subsample.uniform <- function(InputSample, n.required, n.quantiles = 4, use.QS = T){
	library(dplyr)
	library(Hmisc)

	Hospitals <- InputSample %>%
		transmute(Code, Region, Beds = Hmisc::cut2(Beds, g = n.quantiles), QS = if (use.QS) QS else 1) # Prepare data by discretizing continuous variables like the number of acute beds and by changing QS to a fixed value if not to be used

	Selected.hospitals <- c()
	Hospitals.temp <- data.frame()

	for (i in 1:n.required) { # Until n.required is reached..
		if (nrow(Hospitals.temp) == 0) { # If no hospitals are already selected store all sample hospital here
			Hospitals.temp <- Hospitals %>%
				filter(!(Code %in% Selected.hospitals)) %>% # Remove already selected hospitals
				sample_frac() %>% # Permute order, useful only if QS is not used
				arrange(QS)
		}

		# Extract the first hospital of the temporary list and add it to the list of selected hospitals
		Extracted.hospital <- Hospitals.temp[1,]
		Selected.hospitals <- c(Selected.hospitals, Extracted.hospital$Code)

		# Remove from the temporary list all hospitals in the same location/size block of the extracted hospital
		Hospitals.temp <- Hospitals.temp %>%
			filter(!(Region %in% Extracted.hospital$Region & Beds %in% Extracted.hospital$Beds))
	}

	InputSample %>% filter(Code %in% Selected.hospitals) # Filter the initial data by the selected hospital codes
}


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