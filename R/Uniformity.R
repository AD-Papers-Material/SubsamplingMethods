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

uniform.sampling <- function(Input.Sample, n.required, n.quantiles = 4, use.QS = T){

	library(dplyr)
	library(Hmisc)

	# Prepare data by discretizing continuous variables like the number of acute
	# beds and by changing QS to a fixed value if not to be used
	Hospitals <- Input.Sample %>%
		transmute(
			Code, Region,
			Beds = Hmisc::cut2(Beds, g = n.quantiles),
			QS = if (use.QS) QS else 1)

	Selected.hospitals <- c()
	Candidates <- data.frame()

	for (i in 1:n.required) { # Until n.required is reached..

		# If the candidate list is empty, rebuilt it from the non-selected hospitals
		if (nrow(Candidates) == 0) {
			Candidates <- Hospitals %>%
				# Remove already selected hospitals
				filter(!(Code %in% Selected.hospitals)) %>%
				# Permute order, useful only if QS is not used
				sample_frac() %>%
				# Arrange by QS ascending, best hospitals first
				arrange(QS)
		}

		# Extract the first hospital of the temporary list and add it to the list of
		# selected hospitals
		Extracted.hospital <- Candidates[1,]
		Selected.hospitals <- c(Selected.hospitals, Extracted.hospital$Code)

		# Remove from the temporary list all hospitals in the same location/size block
		# of the extracted hospital
		Candidates <- Candidates %>%
			filter(
				!(Region %in% Extracted.hospital$Region),
				!(Beds %in% Extracted.hospital$Beds)
			)
	}

	# Filter the initial data by the selected hospital codes
	Input.Sample %>% filter(Code %in% Selected.hospitals)
}