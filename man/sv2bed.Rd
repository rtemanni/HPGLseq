\name{sv2bed}
\alias{sv2bed}
\title{convert results from DE table to bed track for visualization in IGV}
\usage{
  sv2bed(res, cond1, cond2, Design.file, lfct = 1)
}
\arguments{
  \item{res}{DE results table}

  \item{Cond1}{}

  \item{Cond2}{}

  \item{Design.file}{name of the design file}

  \item{lfct}{log fold-change threshold}
}
\value{
  Bed file with DE genes results
}
\description{
  convert results from DE table to bed track for
  visualization in IGV
}

