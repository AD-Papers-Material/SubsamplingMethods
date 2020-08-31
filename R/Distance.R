#' Get distance between two regions
#'
#' Compute a measure of distance between two regions based on five macroregions.
#' If two hospitals are in the same region the distance is zero and is one they
#' are in the same macroregion. The distance is between 2 and 5 between
#' different macroregions.
#'
#' @param reg1 A vector of regions.
#' @param reg2 A vector of regions.
#' @param zones A data frame which associate regions and macroregions.
#'
#' @import dplyr
#' @import stringr
#'
#' @return A distance score.
#' @export
#'
#' @examples
#' get.location.distance('Regione Sicilia', 'Regione Piemonte')
#' # returns 4

get.location.distance <- function(reg1, reg2, zones = read.csv(file.path('Data', 'Zones.csv'))) {

	library(dplyr)
	library(stringr)

	if (length(reg1) != length(reg2)) stop('Region lists should have equal length!')

	# Normalize region name casing
	reg1 <- str_to_lower(reg1)
	reg2 <- str_to_lower(reg2)
	zones$Region <- zones$Region %>% str_to_lower()

	# Pairwise comparisons
	lapply(1:length(reg1), function(i){
		row <- c(reg1[i], reg2[i])

		# Zero if same region
		if (row[1] == row[2]) result <- 0
		# Special case for the distance procedure
		else if ('unassigned' %in% row) result <- 1
		else {
			zones <- c(zones$Zone[zones$Region == row[1]], zones$Zone[zones$Region == row[2]])

			# One if same zone
			if (zones[1] == zones[2]) result <- 1

			# The following indicate ad-hoc zone distances
			else if (all(zones %in% c('isole', 'sud'))) result <- 2
			else if (all(zones %in% c('isole', 'centro'))) result <- 3
			else if (all(zones %in% c('isole', 'nord.ovest'))) result <- 4
			else if (all(zones %in% c('isole', 'nord.est'))) result <- 5
			else if (all(zones %in% c('sud', 'centro'))) result <- 2
			else if (all(zones %in% c('sud', 'nord.ovest'))) result <- 3
			else if (all(zones %in% c('sud', 'nord.est'))) result <- 3
			else if (all(zones %in% c('centro', 'nord.ovest'))) result <- 2
			else if (all(zones %in% c('centro', 'nord.est'))) result <- 2
			else if (all(zones %in% c('nord.ovest', 'nord.est'))) result <- 2
		}

		result
	}) %>% unlist
}

#' Compute distance between elements of a one dimensional vector.
#'
#' @param reg1 An ordered factor variable or integer.
#' @param reg2 An ordered factor variable or integer.
#'
#' @import dplyr
#'
#' @return A distance score.
#' @export
#'
#' @examples
#'
#' small <- factor('small', levels = c('small', 'medium', 'big'))
#' big <- factor('big', levels = c('small', 'medium', 'big'))
#' get.onedim.distance(small, big)

get.onedim.distance <- function(reg1, reg2) {

	library(dplyr)

	if (length(reg1) != length(reg2)) stop('Region lists should have equal length!')

	case_when(
		# special case for unassigned hospitals
		as.character(reg1) == 'unassigned' | as.character(reg2) == 'unassigned' ~ 1,
		# simple absolute difference for the rest
		T ~ abs(as.numeric(reg1) - as.numeric(reg2))
	)
}

#' Distance sub-sampling procedure
#'
#' @param Input.Sample A data frame with the initial sample that needs subs-
#'   sampling. It needs to have the columns: Code, Region, Beds, QS.
#' @param Reference.Data Target population level reference data, with
#'   information on number of Beds and Region of all hospitals.
#' @param n.required  The size of the final sub-sample.
#' @param n.quantiles Number of quantiles into which categorize continuous
#'   variables, like the hospital size.
#' @param priority Whether to prioritize adjustment of hospital size or location
#'   bias.
#' @param method Statistical method to evaluate the fit between the subsample
#'   and the target population.
#' @param reallocate Option to perform the random reallocation after the first
#'   subsample is created.
#' @param realloc.steps How many reallocation steps to try.
#' @param steps.to.try Number of steps after which to stop if no useful swaps
#'   have been found. Defaults to 2 * the required subsample size.
#'
#' @import dplyr
#' @import Hmisc
#' @import pbapply
#' @import parallel
#'
#' @return A subset of the data frame passed in Input.Sample with n.required
#'   rows
#' @export
#'
#' @examples
#'
#' # Create a subsample of 55 hospitals prioritizing hospital size, without
#' # reallocation after the initial sampling
#' subsample.distance(Sample.Data, Reference.Data, n.required = 55,
#'     reallocate = F)
#'
#' # This time prioritize location and perform reallocation (Warning: it takes
#' # time!)
#' ## Not run:
#' subsample.distance(Sample.Data, Reference.Data, n.required = 55,
#'     priority = 'location')
#' ## End(Not run)

subsample.distance <- function(Input.Sample, Reference.Data, n.required,
															 n.quantiles = 10, priority = c('size', 'location'),
															 method = c('logLik', 'spearman'), reallocate = T,
															 realloc.steps = 2000, steps.to.try = 2 * n.required){

	library(dplyr)
	library(Hmisc)
	library(pbapply)
	library(parallel)

	priority <- match.arg(priority)
	method <- match.arg(method)

	# Definition of quantiles in the distribution of number of beds according to
	# reference data
	quantiles <- quantile(Reference.Data$Beds, seq(0, 1, length.out = n.quantiles)) %>%
		round()

	Hospitals <- Input.Sample %>%
		transmute(
			Code, QS,

			# Hospital size is quantized according to reference data and location/size
			# blocks are defined
			Region,
			Size.class = Hmisc::cut2(Beds, quantiles),
			Block = paste(Size.class, '-', Region),

			# To each hospital is assigned a fictitious region and size class, the
			# initial value is 'unassigned' which means not selected in the sample
			Region.assigned = 'unassigned',
			Size.class.assigned = factor('unassigned',
																	 levels = c(levels(Size.class), 'unassigned')),

			# The distance associated to the 'unassigned' status is 1
			Distance = 1,
			Size.diff = 1
		)

	# Blocks are identified in the refernce data too
	Reference.Data <- Reference.Data %>% mutate(
		Size.class = Hmisc::cut2(Beds, quantiles),
		Block = paste(Size.class, '-', Region)
	)

	# For each block the expected number or required hospital is computed,
	# translating the target population proportion to the required sample size
	Blocks <- count(Reference.Data, Block, name = 'N.country') %>%
		mutate(
			F.expected = N.country/sum(N.country) * n.required,
			N.expected = round(F.expected),

			# If the required sample size is not reached, additional units are added
			# to blocks closer to be rounded up
			N.expected = if (sum(N.expected) == n.required) N.expected else {
				delta <- (F.expected - N.expected)
				rank <- row_number(-delta) # rank blocks by fractional part
				missing <- 1:(n.required - sum(N.expected))
				case_when(
					rank %in% missing & delta > 0 ~ N.expected + 1,
					T ~ N.expected
				)
			}
		) %>%
		# Join the expected block numerosity with the Reference data, to associate
		# information about region and size class
		left_join(Reference.Data[,c('Region', 'Size.class', 'Block')], by = 'Block') %>%
		distinct %>%
		# Blocks that need more hospital first
		arrange(desc(F.expected))

	## Reservoir of available hospitals to assign
	Available.hospitals <- Hospitals %>%
		slice_sample(prop = 1)

	# First assignment: associate hospitals to block until N_expected is reached;
	# if not enough hospitals for a block are available in the sample, uses
	# hospitals from "similar" blocks.
	Hospitals.reassigned <- pblapply(which(Blocks$N.expected != 0), function(i) {
		Block <- Blocks[i,]

		if (Block$N.expected == 0) return(NULL) # No hospitals required for this block

		# Compute "distances" of the block with all hospitals
		Available.hospitals$Distance <- get.location.distance(
			rep(Block$Region, nrow(Available.hospitals)),
			Available.hospitals$Region
		)
		Available.hospitals$Size.diff <- get.onedim.distance(
			rep(Block$Size.class, nrow(Available.hospitals)),
			Available.hospitals$Size.class
		)

		# Arrange hospitals by distance from the block, prioritizing the
		# chosen characteristic. Hospitals belonging to the block will have
		# zero distance and will be prioritize. QS is used in case of ties
		Selected.hospitals <- Available.hospitals %>% {
			if (priority == 'size') {
				arrange(., Size.diff, Distance, QS)
			} else {
				arrange(., Distance, Size.diff, QS)
			}
		} %>%
			# Select the N_expected hospitals for the block, ordered by distance.
			head(Block$N.expected) %>%
			# "Assign" the block size location to the selected hospitals
			mutate(
				Region.assigned = !!Block$Region,
				Size.class.assigned = !!Block$Size.class
			)

		# Remove the selected hospitals from those still available for sampling
		Available.hospitals <<- Available.hospitals %>%
			filter(!(Code %in% Selected.hospitals$Code))

		Selected.hospitals
	}) %>% bind_rows()

	score <- get.distr.fit(Hospitals.reassigned, Reference.Data, method = method)
	message('First distribution score: ', score)

	# Rebuild the dataset with assigned and not hospitals
	Hospitals <- bind_rows(
		Hospitals.reassigned,
		Hospitals %>% filter(!(Code %in% Hospitals.reassigned$Code))
	)

	if (reallocate) {
		message('Reallocation')

		# Windows doesn't support forking for parallelization, and socketing is
		# actually slower
		options(mc.cores =  if (Sys.info()[['sysname']] == 'Windows') {
			1
		} else parallel::detectCores())

		steps.wo.improvements <- 0

		# Initiate random swapping of hospitals to improve the fit
		pblapply(1:realloc.steps, function(step) {

			steps.wo.improvements <<- steps.wo.improvements + 1

			# iI too many steps are gone without improvement, stops
			if (steps.wo.improvements > steps.to.try) {
				return()
			}

			# Select a random hospital
			Evaluated.hospital <- slice_sample(Hospitals, n = 1)

			# Evaluate swaps with hospital in the opposite assignment status
			Candidate.swaps <- Hospitals %>%
				filter(if (Evaluated.hospital$Region.assigned == 'unassigned') {
					Region.assigned != 'unassigned'
				} else Region.assigned == 'unassigned')

			# Evaluate all candidate swaps, computing the change in distance between
			# the assigned and the actual block
			Best.swap <- mclapply(1:nrow(Candidate.swaps), function(i) {
				Candidates <- bind_rows(Evaluated.hospital, Candidate.swaps[i,])

				# Compute the distances after swapping
				Candidates$swapped.dist <- get.location.distance(
					Candidates$Region, rev(Candidates$Region.assigned)
				)
				Candidates$swapped.size.diff <- get.onedim.distance(
					Candidates$Size.class, rev(Candidates$Size.class.assigned)
				)
				Candidates$swapped.QS <- rev(Candidates$QS)

				data.frame(
					Evaluated = Candidates$Code[1],
					Candidate = Candidates$Code[2],

					# Compute the delta in the distances, a negative value
					# means improvement
					delta.location = with(Candidates, sum(swapped.dist) - sum(Distance)),
					delta.size = with(Candidates, sum(swapped.size.diff) - sum(Size.diff)),
					delta.QS = with(
						filter(Candidates, Region.assigned != 'unassigned'),
						swapped.QS - QS
					)
				)
			}) %>% bind_rows() %>%
				# Arrange the swaps according to priority
				{
					if (priority == 'size') {
						arrange(., delta.size, delta.location, delta.QS)
					} else {
						arrange(., delta.location, delta.size, delta.QS)
					}
				} %>% head(1) # And select the best swap

			# If at least one of the distances real/assigned blocks is imporoved proceed
			if (with(Best.swap, delta.location < 0 | delta.size < 0 | delta.QS < 0)) {

				# Select the hospitals of the best swap and invert their assigned
				# blocks. Then recompute the distances
				Best.swap.hospitals <- Hospitals %>%
					filter(Code %in% unlist(Best.swap[,1:2])) %>%
					mutate(
						across(c(Region.assigned, Size.class.assigned), rev),
						Distance = get.location.distance(Region, Region.assigned),
						Size.diff = get.onedim.distance(Size.class, Size.class.assigned)
					)

				# Remove the swapped hospitals from the sample and add them again with
				# the updated assignment
				Hospital.proposal <- Hospitals %>%
					filter(Code %nin% Best.swap.hospitals$Code) %>%
					rbind(Best.swap.hospitals)

				# Compute the fit with the target population before and after the swap
				score_before <- get.distr.fit(
					filter(Hospitals, Region.assigned != 'unassigned'),
					Reference.Data, method = method)
				score_after <- get.distr.fit(
					filter(Hospital.proposal, Region.assigned != 'unassigned'),
					Reference.Data, method = method)

				# If the the swap didn't decrease the fit, the swap is accepted
				if (score_after >= score_before) {

					message('\nStep: ', step)
					message('Improvement after ', steps.wo.improvements, ' steps')
					message('Score after reassignment:', score_after)

					steps.wo.improvements <<- 0

					Hospitals <<- Hospital.proposal
				}
			}
		})

	}

	# Return the susampled data
	Input.Sample %>%
		filter(Code %in% filter(Hospitals, Region.assigned != 'unassigned')$Code)
}