library(Biostrings)
library(TFBSTools)
library(JASPAR2020)

ID = 'MA0018.4'
output = 'E:/sample/sources/pwm_CREB1.csv'

pfm<-getMatrixByID(JASPAR2020,ID)
pwm<-toPWM(pfm)
write.csv(pwm@profileMatrix, output)

