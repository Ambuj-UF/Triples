
################################################################################################################
# Creates phylogeny triples                                                                                    #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################


library("ape")

extriples <- function (phy, authTaxa, outgroup) {
    treeData <- root(phy, outgroup)
    taxaInc = c()
    for (x in authTaxa) {taxaInc = c(taxaInc, x)}
    taxaInc = c(taxaInc, outgroup)
    if (!FALSE %in% (taxaInc %in% treeData$tip.label)) {
        taxaExc = treeData$tip.label[!treeData$tip.label %in% taxaInc]
        treeData <- drop.tip(treeData, taxaExc)
        
        return(treeData)
    }
        
    else {return("NONE")}
}


triples <- function (input = "", auth = "", outgroup = "", output="Triples.txt") {
    
    stopifnot(input != "")
    stopifnot(auth != "")
    stopifnot(outgroup != "")
    
    out="td"
    options(warn=-1)
    inData = read.tree(input)
    authTaxa = read.table(file=auth,header=FALSE,sep="\n")
    taxaIncPre=c()

    for (i in 1:length(authTaxa$V1)) {
        taxaName = sapply(authTaxa[i,], as.character)
        taxaIncPre = c(taxaIncPre, taxaName)
    }

    taxaInc = combn(taxaIncPre, 3)
    treeDataTrip = list()

    for (counter in 1:length(inData)) {
        treeDataTrip[[counter]] = list()
    }

    for (counter in 1:length(inData)) {
        cat("Creating triples for tree:", counter, "\n")
        store_i = c()
        for (i in 1:length(taxaInc[1,])){
            retData = extriples(inData[[counter]], taxaInc[,i], outgroup)
            if (retData != "NONE") {
                treeDataTrip[[counter]][[i]] = retData
            }
            else {
                store_i = c(store_i, i)
                treeDataTrip[[counter]][[i]] = "NONE"
            }
        }
    }

    sink(output)

    for (i in 1:length(taxaInc[1,])) {
        if (!i %in% store_i) {
            t1 = read.tree(text=paste('(', outgroup, ',', '(', '(', taxaInc[,i][1], ',', taxaInc[,i][2], ')', ',', taxaInc[,i][3], ')', ')', ';', sep=""))
            t2 = read.tree(text=paste('(', outgroup, ',', '(', '(', taxaInc[,i][1], ',', taxaInc[,i][3], ')', ',', taxaInc[,i][2], ')', ')', ';', sep=""))
            t3 = read.tree(text=paste('(', outgroup, ',', '(', '(', taxaInc[,i][2], ',', taxaInc[,i][3], ')', ',', taxaInc[,i][1], ')', ')', ';', sep=""))
            t4 = read.tree(text=paste('(', outgroup, ',', '(', taxaInc[,i][1], ',', taxaInc[,i][2], ',', taxaInc[,i][3], ')', ')', ';', sep=""))

            count1 = 0; count2 = 0; count3 = 0; count4 = 0
            for (counter in 1:length(inData)) {
                if (treeDataTrip[[counter]][[i]] != "NONE"){
                    treeDataTrip[[counter]][[i]]$edge.length = NULL
                    
                    if (all.equal(treeDataTrip[[counter]][[i]], t1) == TRUE) {
                        count1 = count1 + 1
                    }
                    else if (all.equal(treeDataTrip[[counter]][[i]], t2) == TRUE) {
                        count2 = count2 + 1
                    }
                    else if (all.equal(treeDataTrip[[counter]][[i]], t3) == TRUE) {
                        count3 = count3 + 1
                    }
                    else if (all.equal(treeDataTrip[[counter]][[i]], t4) == TRUE) {
                        count4 = count4 + 1
                    }
                }
            }


            if (out == "td") {
                cat(c(paste(paste('(', outgroup, ',', '(', '(', taxaInc[,i][1], ',', taxaInc[,i][2], ')', ',', taxaInc[,i][3], ')', ')', ';', sep=""),
                paste('(', outgroup, ',', '(', '(', taxaInc[,i][1], ',', taxaInc[,i][3], ')', ',', taxaInc[,i][2], ')', ')', ';', sep=""),
                paste('(', outgroup, ',', '(', '(', taxaInc[,i][2], ',', taxaInc[,i][3], ')', ',', taxaInc[,i][1], ')', ')', ';', sep=""),
                paste('(', outgroup, ',', '(', taxaInc[,i][1], ',', taxaInc[,i][2], ',', taxaInc[,i][3], ')', ')', ';', sep=""), sep="\t"), "\n"))
                cat(c(paste(count1, count2, count3, count4, sep="\t")), "\n")
                cat("\n\n")
            }
            
            else {
                cat(c(paste('(', outgroup, ',', '(', '(', taxaInc[,i][1], ',', taxaInc[,i][2], ')', ',', taxaInc[,i][3], ')', ')', ';', sep=""), count1, '\n'))
                cat(c(paste('(', outgroup, ',', '(', '(', taxaInc[,i][1], ',', taxaInc[,i][3], ')', ',', taxaInc[,i][2], ')', ')', ';', sep=""), count2, '\n'))
                cat(c(paste('(', outgroup, ',', '(', '(', taxaInc[,i][2], ',', taxaInc[,i][3], ')', ',', taxaInc[,i][1], ')', ')', ';', sep=""), count3, '\n'))
                cat(c(paste('(', outgroup, ',', '(', taxaInc[,i][1], ',', taxaInc[,i][2], ',', taxaInc[,i][3], ')', ')', ';', sep="")), count4)
                cat("\n\n")
            }
    
        }
    }

    sink()
}


