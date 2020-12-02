library(seqinr)

# We introduce the two sequence which we need to do the exercises
zika_genome <- read.fasta("zika.fasta")
dengue_genome <- read.fasta("dengue.fasta")

zika <- zika_genome[[1]]
dengue <- dengue_genome[[1]]

###########################################################################################
"
Some genomes have long stretches of either GC-rich or AT-rich sequence. Use
a HMM with two different states (“AT-rich” and “GC-rich”) to infer which state of
the HMM is most likely to have generated each nucleotide position in Zika and
Dengue sequences. For the AT-rich state, set pA= 0.329, pc = 0.301, pG = 0.159,
and pT = 0.211. For the GC-rich state, set pA = 0.181, pC = 0.313, pG = 0.307,
and pT = 0.199. Set the probability of switching from the AT-rich state to the GC-
rich state, or conversely, to be 0.3. Make a plot for each virus in order to see the
change points. Which of both viruses has more change points?
"
###########################################################################################
library(HMM)

create_model <- function(trans_matrix, emis_probs ,init_probs)
{
  "
  We create a model with all the parameters that are required for it.
  "
  possible_obs <- c("a","c","g","t")
  states <- c("AT-rich", "CG- rich")
  model = initHMM(states,possible_obs,init_probs,trans_matrix,emis_probs)
  return(model)
}

plot_distribution_states <- function(model,seq)
{
  "
  We plot the sequence and in the graph we can see which part of the sequence belongs to each state.
  "
  #path <- viterbi(model,observations); path
  path=viterbi(model,seq); print(path)
  path <- ifelse(path == "AT-rich", 1, 0) 
  ts.plot(path)
  return(path)
}

switching_state <- function(path)
{
  "
  
  "
  current_state = path[1]
  numbers_states = 1
  for(position in path[2:length(path)])
  {
    if(current_state != position)
    {
      current_state <- position
      numbers_states = numbers_states + 1
    }
  }
  return(numbers_states)
}
# First row for A-T_rich state
transition_matrix <- matrix(c(0.7,0.3,0.3,0.7),2)
# First row for A-T_rich state. Columns for each possible observation (A,C,G,T)
emision_probs <- matrix(c(0.329,0.301,0.159,0.211,0.181,0.313,0.307,0.199),byrow=TRUE,ncol=4);
# First column for A-T rich state
initial_prob <- c(0.5,0.5)

modelHMM=create_model(transition_matrix, emision_probs, initial_prob)

# ZIKA VIRUS
estimated_model_zika = baumWelch(modelHMM,zika, maxIterations=200)
transition_matrix_estim=matrix(c(0.6277285, 0.3722715,0.3270179, 0.6729821),2)
emision_probs_estim = matrix(c(0.2697527, 0.44759023, 0.05732978, 0.2253272,0.2835506, 0.01736273, 0.49589168, 0.2031950),byrow=TRUE,ncol=4);
zika_HMM2 = create_model(transition_matrix_estim, emision_probs_estim, initial_prob)
path_zika <- plot_distribution_states(zika_HMM2, zika)
switching_zika_state <- switching_state(path_zika); switching_zika_state

# DENGUE VIRUS
estimated_model_dengue = baumWelch(modelHMM,dengue, maxIterations=200)
transition_matrix_estim=matrix(c(0.6252144, 0.3747856,0.3388041, 0.6611959),2)
emision_probs_estim = matrix(c(0.3883114, 0.41478588, 0.03363364, 0.1632691,0.2566152, 0.02232966, 0.46089143, 0.2601637),byrow=TRUE,ncol=4);
dengue_HMM2 = create_model(transition_matrix_estim, emision_probs_estim, initial_prob)
path_dengue <- plot_distribution_states(dengue_HMM2, dengue)
switching_dengue_state <- switching_state(path_dengue); switching_dengue_state

#Which have more switching between states
if(switching_zika_state == switching_dengue_state){print("Both have the same numbers of switching between states")}
if(switching_zika_state > switching_dengue_state) {print("Zika has more numbers of switching between states")
}else {print("Dengue has more numbers of switching between states")}



#Length of the sequences
nzika <- length(zika)
print('Length of Zika virus:'); nzika
nden <- length(dengue)
print('Length of Dengue virus: '); nden

#Vectors for the GC content of the chunks of 100 nucleotides of the viral genomes
gczika <- c()
gcden <- c()
cctzika <- c()
cctden <- c()
contzika <- 1
contden <- 1
pos <- 1

###GC content calculation for each chunk and presence of cct trinucleotide FOR ZIKA VIRUS
while (contzika < nzika) {
  #Creates the chunks of length 100
  if ((contzika+100) < nzika){
    chunk <- zika[contzika:(contzika+99)]#Chunk length 100
    #Caculate GC content and add it to the vector
    gcchunk <- GC(chunk)
    gczika[pos]<-gcchunk
    #Checks presence/absence of cct
    for (i in 1:98){
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctzika[pos] <- 1
        break
      } else {
        cctzika[pos] <- 0
      }
    }
    #continues moving in the while loop
    contzika <- contzika + 100
    pos <- pos+1
    #If it reaches the end and, if the chunk is less than 100, stores it until the end
  } 
  if ((contzika+100) >= nzika){
    chunk <- zika[contzika:nzika]#Chunk of length 100
    #calculate GC content and add it to vector
    gcchunk <- GC(chunk)
    gczika[pos]<-gcchunk
    #Checks presence/absence of cct in the chunk
    fin <- length(chunk)-2
    for (i in 1:fin) {
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctzika[pos] <- 1
        break
      } else {
        cctzika[pos] <- 0
      }
    }
    #continues moving in the while loop
    contzika <- contzika + 100
    pos <- pos+1
  }
}

pos <- 1
###GC content calculation for each chunk and presence of cct trinucleotide FOR DENGUE VIRUS
while (contden < nden) {
  #Creates the chunks of length 100
  if ((contden+100) < nden){
    chunk <- dengue[contden:(contden+99)]#Chunk length 100
    #Caculate GC content and add it to the vector
    gcchunk <- GC(chunk)
    gcden[pos]<-gcchunk
    #Checks presence/absence of cct
    for (i in 1:98){
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctden[pos] <- 1
        break
      } else {
        cctden[pos] <- 0
      }
    }
    #continues moving in the while loop
    contden <- contden + 100
    pos <- pos+1
    #If it reaches the end and, if the chunk is less than 100, stores it until the end
  } 
  if ((contden+100) >= nden){
    chunk <- dengue[contden:nden]#Chunk of length 100
    #calculate GC content and add it to vector
    gcchunk <- GC(chunk)
    gcden[pos]<-gcchunk
    #Checks presence/absence of cct in the chunk
    fin <- length(chunk)-2
    for (i in 1:fin) {
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctden[pos] <- 1
        break
      } else {
        cctden[pos] <- 0
      }
    }
    #continues moving in the while loop
    contden <- contden + 100
    pos <- pos+1
  }
}

####GENERAL LINEAR MODEL - LOGISTIC MODEL (WE ARE RELATING THE VALUES OF A VARIABLE WITH A BINARY RESULT)#######

#PLOT OF THE GC content of each virus and direct comparison of the GC content of each virus along the genome.
plot(gczika, type='l', col='darkred', ylab="GC content (per one)")
lines(gcden, type='l', col='blue4')
legend('bottom', legend=c("Zika virus NC_012532.1", "Dengue virus NC_001477"), col=c("darkred","blue4"), lty=1:1, lwd=2:50, cex=0.8, box.lty=1)

#PLOT OF THE GC content (x axis) vs presence/absence of cct (y axis) for ZIKA AND DENGUE VIRUSES.

#For Zika virus
plot(gczika, cctzika, main="Zika virus", type='p', pch=19, cex=.5, xlab="GC content (per one)", ylab="Presence/Absence of trinucleotide 'cct'")

#For Dengue virus
plot(gcden, cctden, main="Dengue virus", type='p', pch=19, cex=.5, xlab="GC content (per one)", ylab="Presence/Absence of trinucleotide 'cct'")


#This code establishes the coefficients of the general lineal model for the logistic link function for the ZIKA VIRUS
#The linkinf function in this case is g(pi)=b0+b1x1 (sólo hay una variable)
relzika <- glm(cctzika~gczika,family=binomial)
summary(relzika)

#This code establishes the coefficients of the general lineal model for the logistic link function for the ZIKA VIRUS3
#The linkinf function in this case is g(pi)=b0+b1x1 (sólo hay una variable)
relden <- glm(cctden~gcden,family=binomial)
summary(relden)

#After doing this, the values of b0 and b1 for each virus are:
#Zika virus: g(pi)=-1.453(not very significant)+5.802*GCcontent
#Dengue virus: g(pi)=-0.6725+4.2237*GCcontent

#With the coefficients b0 and b1 for each virus, we calculate the probability of finding a cct with a certain GC content in a chunk for different GC contents, and we check the limit at which the relation is significant. 

GCcontents <- c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

#Zika
zikaprob <- c()
for (val in 1:9) {
  prob <- 1/(1+exp(6.537-16.459*GCcontents[val]))
  zikaprob[val] <- prob
}

#Dengue
denprob <- c()
for (val in 1:length(GCcontents)) {
  prob <- 1/(1+exp(3.487-9.078*GCcontents[val]))
  denprob[val] <- prob
}

#Plot the probabilities against the GC content
plot(GCcontents, zikaprob, type='l', col='gold', xlab="GC contents (per one)", ylab="Probability of finding 'cct' for a given GC content")
lines(GCcontents, denprob, type='l', col='darkred', xlab="GC contents (per one)", ylab="Probability of finding 'cct' for a given GC content")
legend('bottom', legend=c("Zika virus NC_012532.1", "Dengue virus NC_001477"), col=c("gold","darkred"), lty=1:1, lwd=2:50, cex=0.8, box.lty=1)


#Then we calculate the probability of finding a ctt with a GC content of 0.5 for each virus(this can be done as a example comparing the values obtained in the upper part). 

#Zika
1/(1+exp(6.537-16.459*0.5))
1/(1+exp(6.537-16.459*0.6))

#Dengue
1/(1+exp(3.487-9.078*0.5))
1/(1+exp(3.487-9.078*0.6))

########################################################################

#Now, the following code will repeat the same steps followed previously, but in order to obtain the probability for the presence/absence of cct in chunks of 150 and 200, in order to see if there is any variation in the probability values, and a clearer relation. 

################################ CHUNKS OF 200 ########################
gczika200 <- c()
gcden200 <- c()
cctzika200 <- c()
cctden200 <- c()
contzika <- 1
contden <- 1
pos <- 1

###GC content calculation for each chunk and presence of cct trinucleotide FOR ZIKA VIRUS
while (contzika < nzika) {
  #Creates the chunks of length 200
  if ((contzika+200) < nzika){
    chunk <- zika[contzika:(contzika+199)]#Chunk length 200
    #Caculate GC content and add it to the vector
    gcchunk <- GC(chunk)
    gczika200[pos]<-gcchunk
    #Checks presence/absence of cct
    for (i in 1:198){
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctzika200[pos] <- 1
        break
      } else {
        cctzika200[pos] <- 0
      }
    }
    #continues moving in the while loop
    contzika <- contzika + 200
    pos <- pos+1
    #If it reaches the end and, if the chunk is less than 100, stores it until the end
  } 
  if ((contzika+200) >= nzika){
    chunk <- zika[contzika:nzika]#Last chunk
    #calculate GC content and add it to vector
    gcchunk <- GC(chunk)
    gczika200[pos]<-gcchunk
    #Checks presence/absence of cct in the chunk
    fin <- length(chunk)-2
    for (i in 1:fin) {
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctzika200[pos] <- 1
        break
      } else {
        cctzika200[pos] <- 0
      }
    }
    #continues moving in the while loop
    contzika <- contzika + 200
    pos <- pos+1
  }
}

pos <- 1
###GC content calculation for each chunk and presence of cct trinucleotide FOR DENGUE VIRUS
while (contden < nden) {
  #Creates the chunks of length 200
  if ((contden+200) < nden){
    chunk <- dengue[contden:(contden+199)]#Chunk length 200
    #Caculate GC content and add it to the vector
    gcchunk <- GC(chunk)
    gcden200[pos]<-gcchunk
    #Checks presence/absence of cct
    for (i in 1:198){
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctden200[pos] <- 1
        break
      } else {
        cctden200[pos] <- 0
      }
    }
    #continues moving in the while loop
    contden <- contden + 200
    pos <- pos+1
    #If it reaches the end and, if the chunk is less than 200, stores it until the end
  } 
  if ((contden+200) >= nden){
    chunk <- dengue[contden:nden]#Last chunk
    #calculate GC content and add it to vector
    gcchunk <- GC(chunk)
    gcden200[pos]<-gcchunk
    #Checks presence/absence of cct in the chunk
    fin <- length(chunk)-2
    for (i in 1:fin) {
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctden200[pos] <- 1
        break
      } else {
        cctden200[pos] <- 0
      }
    }
    #continues moving in the while loop
    contden <- contden + 200
    pos <- pos+1
  }
}

# Now we calculate the parameters beta 0 and 1 of the general linear model for each virus. 
relzika <- glm(cctzika200~gczika200,family=binomial)
summary(relzika)

relden<- glm(cctden200~gcden200,family=binomial)
summary(relden)

#With the coefficients b0 and b1 for each virus, we calculate the probability of finding a cct with a certain GC content in a chunk for different GC contents, and we check the limit at which the relation is significant. 
#Zika
zikaprob200 <- c()
for (val in 1:9) {
  prob <- 1/(1+exp(24.52-60.74*GCcontents[val]))
  zikaprob200[val] <- prob
}

#Dengue
denprob200 <- c()
for (val in 1:length(GCcontents)) {
  prob <- 1/(1+exp(10.32-25.85*GCcontents[val]))
  denprob200[val] <- prob
}

############################################# CHUNKS OF 150 #######################################################################

#Vectors for the GC content of the chunks of 100 nucleotides of the viral genomes
gczika150 <- c()
gcden150 <- c()
cctzika150 <- c()
cctden150 <- c()
contzika <- 1
contden <- 1
pos <- 1

###GC content calculation for each chunk and presence of cct trinucleotide FOR ZIKA VIRUS
while (contzika < nzika) {
  #Creates the chunks of length 150
  if ((contzika+150) < nzika){
    chunk <- zika[contzika:(contzika+149)]#Chunk length 150
    #Caculate GC content and add it to the vector
    gcchunk <- GC(chunk)
    gczika150[pos]<-gcchunk
    #Checks presence/absence of cct
    for (i in 1:148){
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctzika150[pos] <- 1
        break
      } else {
        cctzika150[pos] <- 0
      }
    }
    #continues moving in the while loop
    contzika <- contzika + 150
    pos <- pos+1
    #If it reaches the end and, if the chunk is less than 150, stores it until the end
  } 
  if ((contzika+150) >= nzika){
    chunk <- zika[contzika:nzika]#Last chunk
    #calculate GC content and add it to vector
    gcchunk <- GC(chunk)
    gczika150[pos]<-gcchunk
    #Checks presence/absence of cct in the chunk
    fin <- length(chunk)-2
    for (i in 1:fin) {
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctzika150[pos] <- 1
        break
      } else {
        cctzika150[pos] <- 0
      }
    }
    #continues moving in the while loop
    contzika <- contzika + 150
    pos <- pos+1
  }
}

pos <- 1
###GC content calculation for each chunk and presence of cct trinucleotide FOR DENGUE VIRUS
while (contden < nden) {
  #Creates the chunks of length 150
  if ((contden+150) < nden){
    chunk <- dengue[contden:(contden+149)]#Chunk length 150
    #Caculate GC content and add it to the vector
    gcchunk <- GC(chunk)
    gcden150[pos]<-gcchunk
    #Checks presence/absence of cct
    for (i in 1:148){
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctden150[pos] <- 1
        break
      } else {
        cctden150[pos] <- 0
      }
    }
    #continues moving in the while loop
    contden <- contden + 150
    pos <- pos+1
    #If it reaches the end and, if the chunk is less than 150, stores it until the end
  } 
  if ((contden+150) >= nden){
    chunk <- dengue[contden:nden]#Last chunk
    #calculate GC content and add it to vector
    gcchunk <- GC(chunk)
    gcden150[pos]<-gcchunk
    #Checks presence/absence of cct in the chunk
    fin <- length(chunk)-2
    for (i in 1:fin) {
      if (chunk[i]=="c" & chunk[i+1]=="c" & chunk[i+2]=="t"){
        cctden150[pos] <- 1
        break
      } else {
        cctden150[pos] <- 0
      }
    }
    #continues moving in the while loop
    contden <- contden + 150
    pos <- pos+1
  }
}

# Now we calculate the parameters beta 0 and 1 of the general linear model for each virus. 
relzika <- glm(cctzika150~gczika150,family=binomial)
summary(relzika)

relden<- glm(cctden150~gcden150,family=binomial)
summary(relden)

#With the coefficients b0 and b1 for each virus, we calculate the probability of finding a cct with a certain GC content in a chunk for different GC contents, and we check the limit at which the relation is significant. 
#Zika
zikaprob150 <- c()
for (val in 1:9) {
  prob <- 1/(1+exp(11.437-29.915*GCcontents[val]))
  zikaprob150[val] <- prob
}

#Dengue
denprob150 <- c()
for (val in 1:length(GCcontents)) {
  prob <- 1/(1+exp(1.917-6.964*GCcontents[val]))
  denprob150[val] <- prob
}

#COMPARING THE PROBABILITY OF BEING CCT PRESENT IN A CHUNK OF A PARTICULAR LENGTH WITH EACH CONTENT OF GC
#Zika virus
plot(GCcontents, zikaprob200, main="Zika virus", type='l', col='red', xlab="GC contents (per one)", ylab="Probability of finding 'cct' for a given GC content") #Chunks of 200
lines(GCcontents, zikaprob150, type='l', col='gold') #Chunks of 150
lines(GCcontents, zikaprob, type='l', col='blue')#Chunks of 100
legend('bottom', legend=c("Chunks of 200", "Chunks of 150", "Chunks of 100"), col=c("red","gold", "blue"), lty=1:1, lwd=2:50, cex=0.8, box.lty=1)


#Dengue virus
plot(GCcontents, denprob200, main="Dengue virus", type='l', col='red', xlab="GC contents (per one)", ylab="Probability of finding 'cct' for a given GC content") #Chunks of 200
lines(GCcontents, denprob150, type='l', col='gold') #Chunks of 150
lines(GCcontents, denprob, type='l', col='blue') #Chunks of 100
legend('bottom', legend=c("Chunks of 200", "Chunks of 150", "Chunks of 100"), col=c("red","gold", "blue"), lty=1:1, lwd=2:50, cex=0.8, box.lty=1)


#The two virus together
plot(GCcontents, zikaprob200, main="Zika virus and Dengue virus", type='l', col='red', xlab="GC contents (per one)", ylab="Probability of finding 'cct' for a given GC content") 
lines(GCcontents, denprob200, type='l', col='blue') 
lines(GCcontents, zikaprob150, type='l', col='brown1')
lines(GCcontents, denprob150, type='l', col='cornflowerblue') 
lines(GCcontents, zikaprob, type='l', col='coral')
lines(GCcontents, denprob, type='l', col='cyan')
legend('bottom', legend=c("Zika - chunks of 200", "Zika - chunks of 150", "Zika - chunks of 100", "Dengue - chunks of 200", "Dengue - chunks of 150", "Dengue - chunks of 100"), col=c("red","brown1", "coral", "blue", "cornflowerblue", "cyan"), lty=1:1, lwd=2:50, cex=0.8, box.lty=1)

#SOME PROBABILITES OF 0.5 GC content
1/(1+exp(24.52-60.74*0.5)) #Zika 200
1/(1+exp(10.32-25.85*0.5)) #Dengue 200
1/(1+exp(11.437-29.91*0.5)) #Zika 150
1/(1+exp(1.917-6.964*0.5)) #Dengue 150





















