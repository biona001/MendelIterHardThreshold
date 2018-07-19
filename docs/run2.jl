using IHT, PLINK, RegressionTools
(xbed, ybed) = read_plink_data("test")
kevin_result = L0_reg(xbed, ybed, 2)