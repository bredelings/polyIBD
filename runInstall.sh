#!/usr/bin/env bash

numThreads=1
if [ "$#" -eq 1 ]; then
   numThreads=$1
fi

echo "Rcpp::compileAttributes()"| R --vanilla --slave
echo "roxygen2::roxygenise()"| R --vanilla --slave
#echo "install.packages(\"./\", type = \"source\", repos = NULL, NCPUS=$numThreads)" | R --vanilla --slave
echo "devtools::load_all()" | R --vanilla --slave
