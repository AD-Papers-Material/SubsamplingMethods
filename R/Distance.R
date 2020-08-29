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
		else if ('deposit' %in% row) result <- 1
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
		# special case for the distance procedure
		as.character(reg1) == 'deposit' | as.character(reg2) == 'deposit' ~ 1,
		T ~ abs(as.numeric(reg1) - as.numeric(reg2))
	)
}

subsample.distance <- function(Input.Sample, Reference.Data, n.required,
															 n.quantiles = 10, priority = c('size', 'location'),
															 fit.method = c('logLik', 'spearman')){

	library(dplyr)
	library(Hmisc)

	priority <- match.arg(priority)
	method <- match.arg(method)

	quantiles <- quantile(Reference.Data$Beds, seq(0, 1, length.out = n.quantiles)) %>%
		round()

	# split.points <- Reference.Data$Beds %>% Hmisc::cut2(g = n.quantiles, onlycuts = T)
	#
	# class.levels <- cut.by.splits(Beds.Data$Dimensione, split.points, ordered = F) %>% levels() %>% c('deposit')


	message('Starting point')

	Hospitals <- Input.Sample %>%
		transmute(
			Code, Beds, Region,
			Region.assigned = 'deposit',
			Beds = as.numeric.robust(Beds),
			Class.orig = cut.by.splits(Beds, split.points, ordered = F) %>% factor(levels = class.levels),
			Class.actual = 'deposit',
			Class.ref = cut.by.splits(Beds, ecdc.cut.points, ordered = F),
			HAI.risk = if ('HAI.risk' %in% names(.)) HAI.risk else NULL,
			HAI.risk.QS = if ('HAI.risk.QS' %in% names(.)) HAI.risk.QS else NULL,
			Quality.score = as.numeric.robust(Quality.score),
			Distance = 1,
			Size.diff = 1
		)

	Beds.Data <- Beds.Data %>% mutate(
		Class = cut.by.splits(Dimensione, split.points, ordered = F),
		Class.ref = cut.by.splits(Dimensione, ecdc.cut.points, ordered = F),
		HAI.risk = if ('HAI.risk' %in% names(.)) HAI.risk else NULL,
		HAI.risk.QS = if ('HAI.risk.QS' %in% names(.)) HAI.risk.QS else NULL
	)

	Regions <- as.data.frame(table(Beds.Data$Class, Beds.Data$Descrizione.Regione) / nrow(Beds.Data) * 55, stringsAsFactors = F) %>%
		mutate(
			Prob = Freq/sum(Freq),
			N = round(Freq),
			mod = modify_if(Freq %% 1 - .5, ~ .x > 0, ~ NA),
			order = row_number(-mod),
			N = if (sum(N) < to.select) case_when(order %in% 1:(to.select - sum(N)) ~ N + 1, T ~ N) else N
		) %>%
		select(Class = Var1, Region = Var2, Theor.size = N, Prob) %>%
		mutate(Class = factor(Class, levels = class.levels)) %>%
		arrange(desc(Prob))

	# starting.score <- compute.score.summary(Hospitals %>% mutate(Region.actual = Region.orig, Class.actual = Class.orig), Beds.Data, Regions)
	#
	# print(starting.score %>% as.data.frame())

	codes <- Hospitals$Code ## Leave it! will be updated by the loop!

	message('First distribution')
	Hospitals.reassigned <- pblapply(1:nrow(Regions), function(i) {
		region <- Regions[i,]

		if (region$Theor.size == 0) return(NULL)

		hospitals.in.class <- Hospitals %>% filter(Code %in% codes)

		hospitals.in.class$Distance <- get.distance(rep(region$Region, nrow(hospitals.in.class)), hospitals.in.class$Region.orig)
		hospitals.in.class$Size.diff <- get.size.diff(rep(region$Class, nrow(hospitals.in.class)), hospitals.in.class$Class.orig)

		hospitals.in.class <- hospitals.in.class %>%
			#arrange(Distance, Size.diff, Quality.score) %>%
			arrange(Size.diff, Distance, Quality.score) %>%
			head(region$Theor.size) %>%
			mutate(Region.actual = region$Region, Class.actual = region$Class)

		codes <<- setdiff(codes, hospitals.in.class$Code)

		hospitals.in.class
	}) %>% bind_rows()

	Hospitals.reassigned <- bind_rows(Hospitals.reassigned, Hospitals %>% filter(Code %nin% Hospitals.reassigned$Code)) %>%
		mutate_at(vars(Class.orig, Class.actual), ~ factor(., levels = class.levels))

	#print(compute.score.summary(Hospitals.reassigned, Beds.Data, Regions) %>% as.data.frame())

	# Hospitals.unif <- Hospitals %>% mutate(Region = Region.orig, Class = Class.orig, Region.actual = Region.orig, Class.actual = Class.orig) %>% sample.hospitals.uniform(to.select = to.select)
	# Hospitals.prob <- sample.hospitals.probability(Sample.Hospitals, Beds.Data, to.select = to.select) %>% mutate(Class.ref = cut.by.splits(Beds, ecdc.cut.points, ordered = F))

	# rbind(
	# 	get.cor.score(Hospitals.reassigned, Beds.Data, Regions, method = 'all'),
	# 	get.cor.score(Hospitals.unif, Beds.Data, Regions, method = 'all')
	# ) %>%
	# 	do(rbind(., summarise_all(., ~ .[1] - .[2]))) %>% print

	starting.score <- compute.score.summary(Hospitals.reassigned, Beds.Data)
	print(starting.score)

	message('Reallocation')

	steps <- 2000

	steps.wo.improvements <- 0

	pb <- txtProgressBar(min = 0, max = steps, initial = 1, style = 3)

	for (step in 1:steps) {

		steps.wo.improvements <- steps.wo.improvements + 1

		if (steps.wo.improvements > to.select) {
			cat('\n')
			print(step)
			print('Break')
			break()
		}

		first.candidate <- sample_n(Hospitals.reassigned, 1)

		second.candidates <- Hospitals.reassigned %>%
			filter(
				# Region.actual != first.candidate$Region.actual
				# & Class.actual != first.candidate$Class.actual
				if (first.candidate$Region.actual == 'deposit') Region.actual != 'deposit' else Region.actual == 'deposit'
			)

		cases <- lapply(1:nrow(second.candidates), function(i) {
			candidates <- bind_rows(first.candidate, second.candidates[i,])

			candidates$swapped.dist <- get.distance(candidates$Region.orig, rev(candidates$Region.actual))
			candidates$swapped.size.diff <- get.size.diff(candidates$Class.orig, rev(candidates$Class.actual))
			candidates$swapped.quality <- rev(candidates$Quality.score)

			data.frame(
				first = candidates$Code[1],
				second = candidates$Code[2],
				first.actual = candidates$Region.actual[1],
				second.actual = candidates$Region.actual[2],
				first.theoric = candidates$Region.orig[1],
				second.theoric = candidates$Region.orig[2],
				delta.dist = sum(candidates$swapped.dist - candidates$Distance),
				delta.size.dist = sum(candidates$swapped.size.diff - candidates$Size.diff),
				delta.quality = with(candidates[candidates$Region.actual != 'deposit',], sum(swapped.quality) - sum(Quality.score)) # Actually I don't remember why I used a different code for the sum from the one used for the distance
			)
		}) %>% bind_rows()

		#cases <- cases %>% arrange(delta.dist, delta.size.dist, delta.quality) %>% head(1)
		cases <- cases %>% arrange(delta.size.dist, delta.dist, delta.quality) %>% head(1)

		candidates <- Hospitals.reassigned %>% filter(Code %in% c(first.candidate$Code, cases$second)) %>%
			mutate_at(vars(Region.actual, Class.actual, Quality.score), rev) %>%
			mutate(
				Distance = get.distance(Region.orig, rev(Region.actual)),
				Size.diff = get.size.diff(Class.orig, rev(Class.actual))
			)

		score_before <- get.cor.score(Hospitals.reassigned, Beds.Data, method = mtd)

		Hospital.proposal <- Hospitals.reassigned %>% filter(Code %nin% candidates$Code) %>% rbind(candidates)

		score_after <- get.cor.score(Hospital.proposal, Beds.Data, method = mtd)

		condition <- case_when(
			score_after > score_before ~ F,
			score_after < score_before ~ T,
			cases$delta.dist >= 0 & cases$delta.size.dist >= 0 & cases$delta.quality >= 0 ~ F,
			T ~ T
		)

		# score.diff <- rbind(
		# 	get.cor.score(Hospitals.reassigned, Beds.Data, Regions, all = T),
		# 	get.cor.score(Hospital.proposal, Beds.Data, Regions, all = T)
		# ) %>%
		# 	{
		# 		DF <- .
		# 		set_colnames(DF, DF[1,] %>% str_remove(':.*'))
		# 	} %>%
		# 	as.data.frame() %>%
		# 	mutate_all(~ str_remove(., '.* ') %>% as.numeric()) %>%
		# 	summarise_all(~ diff(.))
		#
		# score.diff.vec <- bind_rows(
		# 	if (exists('score.diff.vec')) score.diff.vec else NULL,
		# 	score.diff
		# )

		# score.vec <- bind_rows(
		# 	if (exists('score.vec')) score.vec else NULL,
		# 	get.cor.score(Hospital.proposal, Beds.Data, Regions, all = T) %>% lapply(function(x) {
		# 		nm <- str_remove(x, ':.*')
		# 		val <- str_remove(x, '.* ') %>% as.numeric()
		# 		data.frame(val) %>% set_colnames(nm)
		# 	}) %>% bind_cols()
		# )

		if (condition) {

			#print(steps.wo.improvements)

			steps.wo.improvements <- 0

			Hospitals.reassigned <- Hospital.proposal

			# Hospitals[c(first.candidate$Code, cases$second),]$Region.actual <- rev(candidates$Region.actual)
			# Hospitals[c(first.candidate$Code, cases$second),]$Distance <- get.distance(rev(candidates$Region.actual), candidates$Region.orig)
		}

		setTxtProgressBar(pb, step)
	}

	Hospitals.reassigned <- Hospitals.reassigned %>% filter(Region.actual != 'deposit')

	# res <- rbind(
	# 	get.cor.score(Hospitals.reassigned, Beds.Data, Regions, method = 'all'),
	# 	get.cor.score(Hospitals %>% mutate(Region.actual = Region.orig, Class.actual = Class.orig), Beds.Data, Regions, method = 'all')
	# ) %>%
	# 	summarise_all(., ~ .[1] / .[2]) %>%
	# 	mutate(
	# 		pop.risk = mean(Beds.Data$HAI.risk), pop.risk.w = with(Beds.Data, weighted.mean(HAI.risk, Dimensione)),
	# 		sample.risk = mean(Sample.Hospitals$HAI.risk), sample.risk.w = with(Sample.Hospitals, weighted.mean(HAI.risk, Beds)),
	# 		subsample.risk = mean(Hospitals.reassigned$HAI.risk), subsample.risk.w = with(Hospitals.reassigned, weighted.mean(HAI.risk, Beds)),
	# 		method = method, groups = n_groups, steps = step, score_change = 1 - score_after/score_before
	# 	)
	#
	# print(res)
	#
	# return(res)

	final.score <- compute.score.summary(Hospitals.reassigned, Beds.Data)

	message('After reassignment:')
	print(final.score)

	Scores <- rbind(
		final.score %>% mutate(Status = 'Final'),
		starting.score %>% mutate(Status = 'Starting')
	)

	Scores <- Scores %>%
		group_by(Class.ref) %>%
		summarise_if(is.numeric, ~ sprintf('%g -> %g (%s)', .[2], .[1], percent(.[1]/.[2] - 1))) %>%
		lapply(function(x) x %>% set_names(Scores$Class.ref %>% levels)) %>%
		tail(-1)
	# rbind(
	# 	get.cor.score(Hospitals.reassigned, Beds.Data, Regions, method = 'all'),
	# 	get.cor.score(Hospitals.rnd, Beds.Data, Regions, method = 'all'),
	# ) %>%
	# 	do(rbind(., summarise_all(., ~ .[1] - .[2])))

	# message(sprintf('Distance: %d (%.3g)', sum(Hospitals.final$Distance[Hospitals.final$Region.actual != 'deposit']), sum(Hospitals.final$Distance[Hospitals.final$Region.actual != 'deposit']) / sum(Hospitals$Distance[Hospitals$Region.actual != 'deposit'])))
	# message(sprintf('Quality: %.4g (%.3g)', sum(Hospitals.final$Quality.score[Hospitals.final$Region.actual != 'deposit']), sum(Hospitals.final$Quality.score[Hospitals.final$Region.actual != 'deposit']) / sum(Hospitals$Quality.score[Hospitals$Region.actual != 'deposit'])))
	# message('Assigned: ', Hospitals %>% filter(Region.actual != 'deposit') %>% nrow)
	# message(sprintf('Travelers: %d (%.3g)', Hospitals.final %>% filter(Region.actual != 'deposit' & Region.actual != Region.orig) %>% nrow, nrow(Hospitals.final %>% filter(Region.actual != 'deposit' & Region.actual != Region.orig)) / nrow(Hospitals %>% filter(Region.actual != 'deposit' & Region.actual != Region.orig))))
	# print(
	# 	(ggplot(Hospitals, aes(Region.orig)) + ggtitle('Sample')) /
	# 		(ggplot(Beds.Data, aes(Descrizione.Regione)) + ggtitle('Ref')) /
	# 		(ggplot(Hospitals.reassigned, aes(Region.orig)) + ggtitle('Sample')) &
	# 		geom_bar(aes(y = stat(count) / sum(stat(count)), fill = Class.ref), position = position_dodge(preserve = 'single'), show.legend = F) &
	# 		theme.obligue.x.axis
	# )
	#
	# data.frame(
	# 	Sample = with(Hospitals, table(Region.orig, Class.ref)) %>% as.data.frame(),
	# 	Reg = with(Beds.Data, table(Descrizione.Regione, Class.ref)) %>% as.data.frame() %>% select(3),
	# 	Subsample = with(Hospitals.reassigned, table(Region.orig, Class.ref)) %>% as.data.frame() %>% select(3)
	# ) %>% mutate_at(vars(matches('Freq'), ~ . / sum(.))) %>% View
	#
	#
	list(
		Hospitals = Hospitals.reassigned,
		Scores = Scores
	)
}