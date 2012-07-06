fastq<- read.table('/Users/liliahedberg/Documents/U2OS/Allele_specific_expression/SNP_list_U2OS_allelespec.unix.fastq')

#fastq<- read.table('/Users/liliahedberg/Desktop/release.0.0/simulate.reads/ouput.example.fastq')


#You need to count the rows in your fastq file and change the iterator.
temp<-cbind(as.character(fastq[seq(1,1904397,4),1]), as.character(fastq[seq(2,1904398,4),1]), as.character(fastq[seq(3,1904399,4),1]), as.character(fastq[seq(4,1904400,4),1]))

make.mutations <- function(temp, average.rate) {
	
	muted<-apply(temp, 1, function(i){
		string <- i[2]
		vect <- strsplit(string, split="")
		
		#print(vect)
		nmut<- rbinom(1,74, rbeta(1, 100*average.rate, 100-(100*average.rate))) #The second parameter in rbinom should be Readlength-1 
		
		if(nmut > 0){
			
			for(j in 1:nmut){
				mut.base <- floor(runif(1,1,75.999)) #The third parameter in runif should be readlength + 0.999
				#print(mut.base)
				#print(vect[[1]][mut.base])
				
				if(toupper(vect[[1]][mut.base]) == "A"){
					vect[[1]][mut.base] <- sample(c("C","G", "T"), 1)
				} else
					
					if(toupper(vect[[1]][mut.base]) == "C"){
						vect[[1]][mut.base] <- sample(c("A","G", "T"), 1)
					} else
						
						if(toupper(vect[[1]][mut.base]) == "G"){
							vect[[1]][mut.base] <- sample(c("C","A", "T"), 1)
						} else
							if(toupper(vect[[1]][mut.base]) == "T"){
								vect[[1]][mut.base] <- sample(c("C","G", "A"), 1)
							}
						}
					}
					
			return(paste(vect[[1]][1],vect[[1]][2],vect[[1]][3],vect[[1]][4],vect[[1]][5],vect[[1]][6],vect[[1]][7],vect[[1]][8],vect[[1]][9],vect[[1]][10],vect[[1]][11],vect[[1]][12],vect[[1]][13],vect[[1]][14],vect[[1]][15],vect[[1]][16],vect[[1]][17],vect[[1]][18],vect[[1]][19],vect[[1]][20],vect[[1]][21],vect[[1]][22],vect[[1]][23],vect[[1]][24],vect[[1]][25],vect[[1]][26],vect[[1]][27],vect[[1]][28],vect[[1]][29],vect[[1]][30],vect[[1]][31],vect[[1]][32],vect[[1]][33],vect[[1]][34],vect[[1]][35],vect[[1]][36],vect[[1]][37],vect[[1]][38],vect[[1]][39],vect[[1]][40],vect[[1]][41],vect[[1]][42],vect[[1]][43],vect[[1]][44],vect[[1]][45],vect[[1]][46],vect[[1]][47],vect[[1]][48],vect[[1]][49],vect[[1]][50],vect[[1]][51],vect[[1]][52],vect[[1]][53],vect[[1]][54],vect[[1]][55],vect[[1]][56],vect[[1]][57],vect[[1]][58],vect[[1]][59],vect[[1]][60],vect[[1]][61],vect[[1]][62],vect[[1]][63],vect[[1]][64],vect[[1]][65],vect[[1]][66],vect[[1]][67],vect[[1]][68],vect[[1]][69],vect[[1]][70],vect[[1]][71],vect[[1]][72],vect[[1]][73],vect[[1]][74],vect[[1]][75], sep=''))
})   #You need to make sure that this is as long as you reads.
	
	return(muted)
}

mut.01.beta<-matrix(t(cbind(temp[,1], make.mutations(temp, 0.01), temp[,3:4])), ncol=1)
write.table(mut.01.beta, file='/Users/liliahedberg/Documents/U2OS/Allele_specific_expression/U2OS_snp_list.unix_mut1.fastq', quote=F, row.names=F, col.names=F)

#write.table(mut.01.beta, file='/Users/liliahedberg/Desktop/release.0.0/simulate.reads/ouput.example.mut1.fastq', quote=F, row.names=F, col.names=F)

mut.02.beta<-matrix(t(cbind(temp[,1],make.mutations(temp, 0.02), temp[,3:4])), ncol=1)
write.table(mut.01.beta, file='/Users/liliahedberg/Documents/U2OS/Allele_specific_expression/U2OS_snp_list.unix_mut2.fastq', quote=F, row.names=F, col.names=F)

#write.table(mut.01.beta, file='/Users/liliahedberg/Desktop/release.0.0/simulate.reads/ouput.example.mut2.fastq', quote=F, row.names=F, col.names=F)

