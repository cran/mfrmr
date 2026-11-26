\name{load_mfrmr_data}
\alias{load_mfrmr_data}
\title{Load Bundled Synthetic Many-Facet Example Datasets}
\description{
Returns one of the synthetic many-facet rating datasets bundled with \pkg{mfrmr} as a ready-to-use \code{data.frame} for learning, testing, and reproducible examples.
}
\usage{
load_mfrmr_data(name = "example_core", details = FALSE)
}
\arguments{
  \item{name}{Dataset key, or \code{"list"} for the available keys. The default is \code{"example_core"}.}
  \item{details}{With \code{name = "list"}, return an annotated catalogue containing counts, intended uses and design notes.}
}
\value{
A data frame for the selected dataset, a character vector of keys, or a catalogue data frame.
The rating datasets contain person/facet identifiers and observed scores. The
\code{example_gtheory_crossed} and \code{example_gtheory_nested} tables contain
separate \code{Content} and \code{Organization} outcomes. The
\code{example_agreement} table identifies units by Person and Task. The
\code{example_operational_design} table is a planned roster without scores.
Inspect column roles before fitting; these different table types are not
interchangeable inputs to \code{fit_mfrm()}.
}
\seealso{
\code{\link{mfrmr_example_data}}, \code{\link{ej2021_data}}, \code{\link{describe_mfrm_data}}
}
\examples{
load_mfrmr_data("list", details = TRUE)
ratings <- load_mfrmr_data("example_operational")
head(ratings)
table(ratings$Rater, ratings$Score)
}
