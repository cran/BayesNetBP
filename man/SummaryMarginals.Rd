% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F2_SummaryMarginal.R
\name{SummaryMarginals}
\alias{SummaryMarginals}
\title{Summary a continuous marginal distribution}
\usage{
SummaryMarginals(marginals)
}
\arguments{
\item{marginals}{the marginal distributions obtained from \code{\link{Marginals}} function}
}
\value{
a \code{data.frame} object containing information about the marginal distributions for continuous variables.
The marginal distributions of continous variables in a CG-BN model are mixtures of Gaussian distributions.
Therefore, besides the mean and standard deviation, the object has an additional column to specify the number of Gaussian
mixtures.

\describe{
 \item{\code{mean}}{the mean value of a Gaussian distribution.}
 \item{\code{sd}}{the standard deviation of a Gaussian distribution.}
 \item{\code{n}}{the number of Gaussian distributions in the mixture.}
}
}
\description{
This function summary the marginal distributions of continuous variables by outputing the
mean, standard deviation, and number of subpopulations
}
\examples{

data(liver)
tree.init.p <- Initializer(dag=liver$dag, data=liver$data,
                           node.class=liver$node.class,
                           propagate = TRUE)
marg <- Marginals(tree.init.p, c("HDL", "Ppap2a", "Neu1"))
SummaryMarginals(marginals=marg)

}
\references{
Cowell, R. G. (2005). Local propagation in conditional Gaussian Bayesian networks.
Journal of Machine Learning Research, 6(Sep), 1517-1550. \cr
\cr
Yu H, Moharil J, Blair RH (2020). BayesNetBP: An R Package for Probabilistic Reasoning in Bayesian
Networks. Journal of Statistical Software, 94(3), 1-31. <doi:10.18637/jss.v094.i03>.
}
\seealso{
\code{\link{Marginals}}
}
