context("test-transprobs")

test_that("`getTransProbs` works", {
  # setup sim parameters
  f <- 0.5 # this is the proportion of genetic relatedness we are trying to infer
  rho <- 1e-7 # this is the recombination rate that is assumed to be known
  k <- 1 # this is the number of generations that we are trying to infer
  m1 <- 1 # this is the multiplicity of infection for sample1 that we are trying to infer
  m2 <- 1 # this is the multiplicity of infection for sample2 that we are trying to infer
  zmax <- max(m1, m2)-1
  snpdist <- 100 # simulate some positions/genomic coordinates
  m_max <- 6 # some high MOI number that is the max for all samples



  E <- getTransProbs(f, rho, k, zmax)


  # Cpp code to R code
  Evalues <-  dplyr::bind_rows(E["Evalues"])
  Evectors <- dplyr::bind_rows(E["Evectors"])
  Esolve <- dplyr::bind_rows(E["Esolve"])

  replicate(n=length(snpdist), simplify=F, expr =
              replicate(n=(zmax+1), simplify=F, expr =
                                      list(x=matrix(0, nrow = (m_max-1), ncol=(m_max-1)))
              )
  )


  })
