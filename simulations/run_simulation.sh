#!/bin/bash

# accuracy

nohup Rscript simulation_accuracy.R --model logistic --dimension 50 --theta dense > ./log/logistic_dense_accuracy.out &

nohup Rscript simulation_accuracy.R --model logistic --dimension 50 --theta null > ./log/logistic_null_accuracy.out &

nohup Rscript simulation_accuracy.R --model logistic --dimension 50 --theta sparse > ./log/logistic_sparse_accuracy.out &

nohup Rscript simulation_accuracy.R --model lm --dimension 50 --theta dense > ./log/lm_dense_accuracy.out &

nohup Rscript simulation_accuracy.R --model lm --dimension 50 --theta null > ./log/lm_null_accuracy.out &

nohup Rscript simulation_accuracy.R --model lm --dimension 50 --theta sparse > ./log/lm_sparse_accuracy.out &

# coverage

nohup Rscript simulation_coverage.R --model logistic --dimension 50 --sample.size 1000000 --theta dense > ./log/logistic_dense_coverage.out &

nohup Rscript simulation_coverage.R --model logistic --dimension 50 --sample.size 1000000 --theta sparse > ./log/logistic_sparse_coverage.out &

nohup Rscript simulation_coverage.R --model logistic --dimension 50 --sample.size 1000000 --theta null > ./log/logistic_null_coverage.out &

nohup Rscript simulation_coverage.R --model lm --dimension 50 --sample.size 1000000 --theta dense > ./log/lm_dense_coverage.out &

nohup Rscript simulation_coverage.R --model lm --dimension 50 --sample.size 1000000 --theta sparse > ./log/lm_sparse_coverage.out &

nohup Rscript simulation_coverage.R --model lm --dimension 50 --sample.size 1000000 --theta null > ./log/lm_null_coverage.out &