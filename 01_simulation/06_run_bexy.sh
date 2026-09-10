#make master file including all idxstats files
ls ../cami_idxstats/* > idxstats.txt 

# to use multiple  threads
#brew install libomp
#export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
#export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"

#run bexy
bexy infer --idxstats cami_idxstats.txt --keepScaffolds  NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,NC_000023.11,NC_000024.10



## visualize results in R

#```{r}
#library(bexy)
#setwd("./BeXY/")

#bex <- bexy("./BeXY/")
#print(bex)
#plot(bex)
#getPosteriorModeSexKaryotypes(bex)

#writePosteriorModeSexKaryotypes(bex, './HMP_bexy_output', threshold_certainty = 0.95)

#summary(bex)

#```

