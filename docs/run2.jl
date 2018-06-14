using IHT, PLINK, RegressionTools
(xbed, ybed) = read_plink_data("test")
output = L0_reg(xbed, ybed, 2)