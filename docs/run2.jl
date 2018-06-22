using IHT, PLINK, RegressionTools
(xbed, ybed) = read_plink_data("gwas 1 data_kevin", delim=',', header=false); output = L0_reg(xbed, ybed, 10)


function run()
    (xbed, ybed) = read_plink_data("gwas 1 data_kevin", delim=',', header=false)
    output = L0_reg(xbed, ybed, 10)
end