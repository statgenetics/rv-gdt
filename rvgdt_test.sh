# rvgdt tests

rvgdt test --geno ./example/rvgdt_test.geno --ped ./example/rvgdt_test.ped --max_iter 100 --analytical

# remove subject
rvgdt test --geno ./example/rvgdt_test.geno --ped ./example/rvgdt_test.ped --max_iter 100 --remove_subject ./example/remove_subjects.txt

# with covariates
rvgdt test --geno ./example/rvgdt_test.geno --ped ./example/rvgdt_test_covariates.ped --max_iter 100 --analytical

# with weights
rvgdt test --geno ./example/rvgdt_test.geno --ped ./example/rvgdt_test_covariates.ped --max_iter 100 --weight ./example/rvgdt_test.weights --analytical