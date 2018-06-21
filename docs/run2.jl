using IHT, PLINK, RegressionTools
(xbed, ybed) = read_plink_data("gwas 1 data_kevin", ',', header=false)
output = L0_reg(xbed, ybed, 10)
