\name{DEAnalysis}
\alias{DEAnalysis}
\title{Core Differential Expression analysis Function}
\usage{
  DEAnalysis(Design.file, counts.dir, DEComparison.file,
    name.prefix, OutputDir, anot.file, source.format,
    SelFeature, fields, adjpval, cleandescription,
    correct.batch, correct.batch.inmodel, filter.lcg)
}
\arguments{
  \item{Design.file}{File with design matrix SampleID /
  Condition / Batch}

  \item{counts.dir}{Directory with counts per gene for each
  sample(one file per sample)}

  \item{DEComparison.file}{File with pairwise comparisons
  to be performed}

  \item{name.prefix}{a name to be appended to generated
  files}

  \item{OutputDir}{output directory}

  \item{anot.file}{annotation file full path}

  \item{source.format}{annotation format G}

  \item{SelFeature}{}

  \item{fields}{}

  \item{adjpval}{}

  \item{cleandescription}{T/F (default= T) change ascii
  code in decrition to char (comon in gff/gtf from
  http://tritrypdb.org/)}

  \item{correct.batch}{}

  \item{correct.batch.inmodel}{}

  \item{filter.lcg}{filter low count genes}
}
\value{
  PCA/ HeatMap /Non Zero GEnes

  annotated DE genes results table

  BED Track with DE genes
}
\description{
  Core Differential Expression analysis Function
}

