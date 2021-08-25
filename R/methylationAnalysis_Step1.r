#' Reads Bionano Cmap files
#'
#' @param cmap  character. Path to Bionano cmap file.
#' @return Data Frame Contains the cmap file information.
#' @examples
#' cmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
#' modcmap <- readingCmap(cmap)
#' @import utils
#' @export


readingCmap <- function(cmap){
        #Reading contig, molecule or reference cmap
        cmap_temp<-file(cmap, "r")
        cmap_temp1<-readLines(cmap_temp, n = -1)
        close(cmap_temp)
        g1 <- grep("#h", cmap_temp1)
        g2 <- grep("#f ",cmap_temp1)
        ##Extracting data
        if (g1+1 == g2){
		    ##reading the cmap file
            dat <- gsub("#h ", "", cmap_temp1)
            #dat4 <- textConnection(dat[g1:length(dat)])
            dat4 <- textConnection(dat[g1])
            r1 <- read.table(dat4, sep = "\t",header=FALSE)
            #dat4 <- textConnection(dat[g1:length(dat)])
            g3<-g1+2
            dat5 <- textConnection(dat[g3:length(dat)])
            r2 <- read.table(dat5, sep = "\t",header=FALSE)
            close(dat5)
            r8 <- data.frame(r2)
            nam <- tryCatch(
                names(r8) <- as.character(unlist(r1)),
                warning = function(w) {
                    print(paste("File:", cmap, "is empty"))
                    return (NA)
                    },
                    error = function(e) {
                    print(paste("File:", cmap, "is empty"))
                    return (NA)
                    }
                )
			###Returning the dara
            if (length(nam) == 0){
                return(NA)
            }
            else{
                names(r8) <- nam
                return(r8)
            }
}
return(r8)
}

#' Reads Bionano xmap files
#'
#' @param xmap character. Path to Bionano xmap file.
#' @return Data Frame Contains the xmap file information.
#' @examples
#' xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
#' modxmap <- readingXmap(xmap)
#' @import utils
#' @export

readingXmap <- function(xmap) {
    ##Reading contig-molecule or reference-contig xmap
    
    xmap = xmap
    conXmap <- file(xmap, "r")
    r10 <- readLines(conXmap, n = -1)
    close(conXmap)
    datfinal <- data.frame()
    g1 <- grep("#h ", r10)
    g2 <- grep("#f ", r10)
	
	##Reading contig xmap
    if (g1+1 == g2){
	   
        dat <- gsub("#h ", "", r10)
        #dat4 <- textConnection(dat[g1:length(dat)])
        dat4 <- textConnection(dat[g1])
        r1 <- read.table(dat4, sep = "\t",header=FALSE)
        #dat4 <- textConnection(dat[g1:length(dat)])
        g3<-g1+2
        dat5 <- textConnection(dat[g3:length(dat)])
        r2 <- read.table(dat5, sep = "\t",header=FALSE)
        #r8<-data.frame(r2)
		nam <- tryCatch(
                names(r2) <- as.character(unlist(r1)),
                warning = function(w) {
                    print(paste("File:", xmap, "is empty"))
                    return (NA)
                    },
                    error = function(e) {
                    print(paste("File:", xmap, "is empty"))
                    return (NA)
                    }
                )
            ##Returning the data
         	if (length(nam) == 0){
                return(NA)
            }
            else{
                names(r2) <- nam
                return(r2)
            }
        
    }
return(r2)
}

#' Identifies and documents the methylation and nick labels
#'
#' @param molcmap character. Path to Bionano molecule cmap.
#' @param xmapdata character. Path to Bionano molecule contig xmap.
#' @param cntNick integer. Nick label value.
#' @param cntMeth integer. Meth label value.
#' @return Data Frame with the molecule cmap, with the exact site number 
#' for methylation and nick labels.
#' @examples
#' xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
#' cmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
#' modcmap <- readingCmap(cmap)
#' modxmap <- readingXmap(xmap)
#' modMolcmap <- modmolcmap(molcmap = modcmap,
#'   xmapdata = modxmap,
#'   cntNick = 1,
#'   cntMeth=2)
#' @import utils
#' @export
modmolcmap <- function(molcmap, xmapdata, cntNick , cntMeth ){
    print("Adding the label location to molecule files")
    molCmapID <- as.numeric(unique(molcmap$CMapId))
    dataContig<-data.frame()
    datf<-data.frame()
    datFinalmol <- data.frame()
    xmapdata = xmapdata
    molcmap = molcmap
    cntNick = cntNick
    cntMeth = cntMeth
    ###counting nick and methylation region in the molcmap
    for (j in seq_along(molCmapID)){
        #print(molCmapID[j])
        dat3<-molcmap[which(molcmap$CMapId==molCmapID[j] & 
            (molcmap$LabelChannel == cntNick 
            | molcmap$LabelChannel == cntMeth)),]
        datnick <- dat3[which(dat3$LabelChannel == cntNick),]
        datmeth <- dat3[which(dat3$LabelChannel == cntMeth),]
        orientMol <- c()
        countnick <- c()
        countmeth <- c()
        molxmap <- xmapdata[which(xmapdata$QryContigID == molCmapID[j]),]
        dat4 <- molcmap[which(molcmap$CMapId==molCmapID[j] 
		    & molcmap$LabelChannel == 0),]
        dat4 <- cbind(dat4[,1:4],
            NickLoc = "-",
            MethLoc = "-",
            OrientationMolecule = molxmap$Orientation,
            dat4[,5:ncol(dat4)],
			row.names = NULL)
        cntnickneg = nrow(datnick) + 1; cntmethneg = nrow(datmeth) + 1
        orient <- unique(as.character(molxmap$Orientation))
        cntnickpos = 0; cntmethpos = 0
            for (hh in 1:nrow(dat3)){
                    if(as.numeric(dat3$LabelChannel[hh])== cntNick){
                cntnickpos = cntnickpos + 1
                #print(paste("hh:", hh, "cntnickpos:", cntnickpos))
                countnick<-c(countnick,cntnickpos)
                orientMol <-  c(orientMol, orient)
                countmeth<-c(countmeth,0)
            }
            else if(as.numeric(dat3$LabelChannel[hh])== cntMeth){
                cntmethpos = cntmethpos + 1
                countmeth<-c(countmeth,cntmethpos)
                orientMol <-  c(orientMol, orient)
                countnick<-c(countnick,0)
            }
            else {
                countmeth<-c(countmeth,0)
                countnick<-c(countnick,0)
            }
        
            
            }
              datt<-cbind(dat3[,1:4],
            NickLoc = countnick,
            MethLoc = countmeth,
            OrientationMolecule = orientMol,
            dat3[,5:ncol(dat3)],
			row.names = NULL)
            datt <- rbind(datt, dat4,row.names = NULL)
        
        datf<-rbind(datf,datt,row.names = NULL)
    }
        'qrycmapData <- cmap
        qrycontig <- as.numeric(unique(cmap$CMapId))
        datf<-data.frame()
        for (j in 1:length(qrycontig)){
            dat3<-qrycmapData[which(qrycmapData$CMapId == qrycontig[j]),]
            cntnick=0;cntmeth=0;countnick<-c();countmeth<-c()
            for (hh in 1:nrow(dat3)){
            if(as.numeric(dat3$LabelChannel[hh])==2){
                cntnick=cntnick+1
        countnick<-c(countnick,cntnick)
        countmeth<-c(countmeth,0)
            }else if(as.numeric(dat3$LabelChannel[hh])==1){
                cntmeth=cntmeth+1
        countmeth<-c(countmeth,cntmeth)
        countnick<-c(countnick,0)
            }else {
                countmeth<-c(countmeth,0)
        countnick<-c(countnick,0)
            }
        }
        datt<-cbind(dat3[,1:4],NickLoc=countnick,MethLoc=countmeth,dat3[,5:ncol(dat3)])
        datf<-rbind(datf,datt)
    }'
    return (datf)
}
#' Untarring the directory if needed
#'
#' @param cmaploc_mol character. Path to Bionano compressed file.
#' @return uncompressed folder with all the files for analysis
#' @examples
#' cmaploc_mol <- system.file("extdata", package="methometR")
#' untarredFolder <-untaringdirectoy(cmaploc_mol)
#' @importFrom utils untar
#' @export
untaringdirectoy <- function(cmaploc_mol){
    print("Untaring directories if untarring option chosen")
    untar(file.path(cmaploc_mol,"molecules.tar.gz"), exdir = cmaploc_mol)
    lmolecules <- list.files(paste0(cmaploc_mol, "molecules/output/contigs/exp_refineFinal1/merged_smaps"), pattern = "*_q.cmap", full.names = TRUE)
    return(lmolecules)
}
#' Mapping nick label location from contig to molecule
#'
#' @param nickrefLoc character. contig cmap data.
#' @param molcmap character. molecule cmap data.
#' @param xmapdata character. contig-molecule xmap data.
#' @param contigID character. query contig ID.
#' @param nickrefLoc character. Path to Bionano compressed file.
#' @param cntNick integer. Nick label value.
#' @param cntMeth integer. Meth label value.
#' @return dataframe containing the molecule file, with the nick position
#' @examples
#' refcontigXmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/ContigRef.xmap", package="methometR")
#' refCmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/hg19ref_r.cmap", package="methometR")
#' Contigqcmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/SampContig_q.cmap", package="methometR")
#' refxmapdat <- readingXmap(refcontigXmap)
#' refcmapdat <- readingCmap(refCmap)
#' contigcmapdat <- readingCmap(Contigqcmap)
#' nickRef<-nickReference(refxmap = refxmapdat, refcmap = refcmapdat, 
#'     contigcmap = contigcmapdat, contigID = 6701, 
#'    returnMethod =c("dataFrame"),  chrom = 4)
#' xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
#' molcmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
#' contigcmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContig_r.cmap", package="methometR")
#' contigID = 6701; cntNick = 1; cntMeth=2
#' modxmap <- readingXmap(xmap)
#' cmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
#' modcmap <- readingCmap(cmap) 
#' modMolcmap <- modmolcmap(molcmap = modcmap,
#'   xmapdata = modxmap,
#'   cntNick = 1,
#'   cntMeth=2)
#' contigmodcmap <- readingXmap(contigcmap)
#' modmolcmapNick <-molNickLoc(nickrefLoc = nickRef,
#'    molcmap = modMolcmap,
#'    xmapdata = modxmap, 
#'    contigID = contigID, cntNick = cntNick, 
#'    cntMeth = cntMeth)
#' @import utils
#' @export
molNickLoc <- function(nickrefLoc,molcmap,xmapdata, contigID, cntNick , 
    cntMeth){
	print("Searching for the nick labels in the molecules")
    'cl <- makeCluster(5) #not to overload your computer
	print(cl)
    registerDoParallel(cl)'
	molcmap <- molcmap
    datf <- molcmap
      xmapdata <- xmapdata
      datFinal <- nickrefLoc
      contigID = contigID
	  cntNick = cntNick
	  cntMeth = cntMeth
	'print(paste("CntNick: ", cntNick))
	  print(paste("cntMeth: ", cntMeth))'
      datFinalmol <- data.frame()
	###Identifying the molecule ID
	#molCmapID<-as.numeric(unique(molcmap$CMapId))
	molCmapID<-as.numeric(unique(molcmap$CMapId))
    #datFinalmol <- foreach (ui = 1:length(molCmapID), .combine=rbind) %dopar% {
	###Extracting the site IDs of the query(moleule) and reference(contig)
	for(ui in 1:length(molCmapID)){    
        #print(ui)
        #print(paste("CmapID",molCmapID[ui]))
        refloc2<-c();queryloc2<-c()
        #qrycontig<-as.numeric(dat1$QryContigID)
        dattt<-xmapdata[which(as.numeric(xmapdata$QryContigID) == molCmapID[ui]),]
        alnt<-dattt$Alignment
        refloc3<-c();queryloc2<-c()
        if(length(alnt) > 1){
            for(ii in length(alnt)){
                stt<-strsplit(as.character(alnt[ii]),split="\\)\\(")
                stt1<-as.character(stt[[1]])
                for (ti in 1:length(stt1)){
                    if(length(grep("\\(",stt1[ti]))==1){
                        st2<-strsplit(stt1[ti],split="\\(")
                        st3<-strsplit(st2[[1]][2],split=",")
                        queryloc2<-c(queryloc2,as.numeric(st3[[1]][2]))
                        refloc3<-c(refloc3,as.numeric(st3[[1]][1]))
                    }
                    else if(length(grep("\\)",stt1[ti]))==1){
                        st2<-strsplit(stt1[ti],split="\\)")
                        st3<-strsplit(st2[[1]][1],split=",")
                        queryloc2<-c(queryloc2,as.numeric(st3[[1]][2]))
                        refloc3<-c(refloc3,as.numeric(st3[[1]][1]))
                    }
                    else{
                        
                        st3<-strsplit(stt1[ti],split=",")
                        queryloc2<-c(queryloc2,as.numeric(st3[[1]][2]))
                        refloc3<-c(refloc3,as.numeric(st3[[1]][1]))
                    }
                }
            }
        }else{
            stt<-strsplit(as.character(alnt),split="\\)\\(")
            stt1<-as.character(stt[[1]])
            for (ti in 1:length(stt1)){
                if(length(grep("\\(",stt1[ti]))==1){
                    st2<-strsplit(stt1[ti],split="\\(")
                    st3<-strsplit(st2[[1]][2],split=",")
                    queryloc2<-c(queryloc2,as.numeric(st3[[1]][2]))
                    refloc3<-c(refloc3,as.numeric(st3[[1]][1]))
                }
                else if(length(grep("\\)",stt1[ti]))==1){
                    st2<-strsplit(stt1[ti],split="\\)")
                    st3<-strsplit(st2[[1]][1],split=",")
                    queryloc2<-c(queryloc2,as.numeric(st3[[1]][2]))
                    refloc3<-c(refloc3,as.numeric(st3[[1]][1]))
                }
                else{
                    
                    st3<-strsplit(stt1[ti],split=",")
                    queryloc2<-c(queryloc2,as.numeric(st3[[1]][2]))
                    refloc3<-c(refloc3,as.numeric(st3[[1]][1]))
                }
            }       
        }
        dattempcont <- data.frame(QueryID = queryloc2, RefID = refloc3)
        ###Extracting based on the refl
        LocRef<-c();chrompos<-c();chromosome<-c();contigsiteID <- c()
        molsiteID <- c(); contig <- c()
            'datttt <- datf[which(datf$CMapId == molCmapID[ui] &
                datf$LabelChannel == cntNick),]
            queryloc3 <- datttt$NickLoc
            unqueryloc <- unique(queryloc3)'
            unqueryloc <- unique(queryloc2)
        for (ll in 1:length(unqueryloc)){
            
            #k1<-as.numeric(hash::keys(h))
            #print(as.numeric(refloc2[ll]))
            #val1<-as.numeric(hash::values(h,keys=as.numeric(refloc2[ll])))
            datp <- dattempcont[which(dattempcont$QueryID == unqueryloc[ll]),]
            val <- c()
            val <- datp$RefID
            #print(val)
            if(length(val)> 1){
                for( ql in 1:length(val)){
                #val <- hash::values(ha6, keys = queryloc3[ll])
                    dattem<-datFinal[which(as.numeric(datFinal$CMapId) == contigID
                        & as.numeric(datFinal$SiteID) == val[ql]),]
                #print(dim(dattem))
                    if(nrow(dattem) >= 1){
                        for(yu in 1:nrow(dattem)){
                            #print(paste("CMapID",cmapID[ii],sep=""))
                            LocRef<-c(LocRef,as.numeric(dattem$Position[yu]))
                            chrompos <- c(chrompos,as.numeric(dattem$chromosomePosition[yu]))
                            contig <- c(contig, contigID)
                            contigsiteID <- c(contigsiteID, dattem$SiteID[yu])
                            molsiteID <- c(molsiteID, unqueryloc[ll])
                            chromosome<-c(chromosome,as.numeric(dattem$chromosome[yu]))
                        }
                    }else{
                        #print(paste("CMapID",cmapID[ii],sep=""))
                        LocRef<-c(LocRef,0)
                        chrompos <- c(chrompos,NA)
                        contig <- c(contig, contigID)
                        contigsiteID <- c(contigsiteID, 0)
                        molsiteID <- c(molsiteID, unqueryloc[ll])
                        chromosome<-c(chromosome,as.numeric(dattem$chromosome[yu]))
                    }
                }
            }else{
                #val <- hash::values(ha6, keys = queryloc3[ll])
                dattem<-datFinal[which(as.numeric(datFinal$CMapId) == contigID
                    & as.numeric(datFinal$SiteID) == val),]
                #print(dim(dattem))
                if(nrow(dattem) >= 1){
                    for(yu in 1:nrow(dattem)){
                        #print(paste("CMapID",cmapID[ii],sep=""))
                        LocRef<-c(LocRef,as.numeric(dattem$Position[yu]))
                        chrompos <- c(chrompos,as.numeric(dattem$chromosomePosition[yu]))
                        contig <- c(contig, contigID)
                        contigsiteID <- c(contigsiteID, dattem$SiteID[yu])
                        molsiteID <- c(molsiteID, unqueryloc[ll])
                        chromosome<-c(chromosome,as.numeric(dattem$chromosome[yu]))
                    }
                }else{
                    #print(paste("CMapID",cmapID[ii],sep=""))
                    LocRef<-c(LocRef,0)
                    chrompos <- c(chrompos,NA)
                    contig <- c(contig, contigID)
                    contigsiteID <- c(contigsiteID, 0)
                    molsiteID <- c(molsiteID, unqueryloc[ll])
                    chromosome<-c(chromosome,chromosome,as.numeric(dattem$chromosome[yu]))
                }
                
            }
        }
        datmol<-data.frame(
            cmapid = rep(as.character(molCmapID[ui]),length(LocRef)),
            chromosome = chromosome,
            chrompos = chrompos,
            LocRef=LocRef, 
            ContigID = contig,
            ContigSiteID = contigsiteID,
            MolsiteID = molsiteID)    
        cmapID<-unique(datmol$cmapid)
        #print(cmapID)
        dattemp1<-datmol[which((datmol$cmapid) == cmapID),]
        datftemp<-datf[which(as.numeric(datf$CMapId) == cmapID),]
        chromo = unique(datmol$chromosome)
        'contigSiteID <- as.numeric(dattemp1$ContigSiteID)
        Position <- (as.numeric(dattemp1$LocRef))'
        molsiteID <- as.numeric(unique(dattemp1$MolsiteID))
        ik = 1; ij = 0
        OrientationMolecule <- as.character(unique(datftemp$OrientationMolecule))
        #pnc <- paste("^", datftemp$NickLoc, "$", sep ="")
        if(OrientationMolecule == "-"){
            len <- nrow(datftemp)
            datFinalmol_rep <- data.frame()
        for(qi in nrow(datftemp):1){
            #print(qi)
            if(datftemp$LabelChannel[qi] == cntNick){
                
                #dattemp3 <- datftemp[which(as.numeric((datftemp$NickLoc)[qi]) == molsiteID[ik]),]
                #print(paste("molsiteID[ik]:",molsiteID[ik]))
                #print(paste("ik",ik))
                if(is.na(molsiteID[ik])) {
                    molsiteID[ik] = 0
                } else{
                    molsiteID[ik] = molsiteID[ik]
                }
                if(datftemp$NickLoc[qi] == as.character(molsiteID[ik])){
                    dattemp2 <- dattemp1[which((dattemp1$MolsiteID) == molsiteID[ik]),]
                    #print(names(dattemp2))
                    contigSiteID <- as.numeric(dattemp2$ContigSiteID)
                    Position <- (as.numeric(dattemp2$LocRef))
                    Contig <- dattemp2$ContigID
                    #print(Contig)
                    chromPosition <- as.numeric(dattemp2$chrompos)
                    chromosome <- as.numeric(dattemp2$chromosome)
                    if(length(contigSiteID)> 1){
					    'print(qi)
						print(OrientationMolecule)'
                        for(lo in 1:length(contigSiteID)){
                            datt<-cbind(datftemp[qi,1:4],
                            chromosome = chromosome[lo],
                            chromPosition = chromPosition[lo],
                            contigPosition = Position[lo], 
                            ContigID = Contig[lo],
                            contigSiteID = contigSiteID[lo], 
                            datftemp[qi,5:ncol(datftemp)],row.names = NULL)
                            datFinalmol_rep<-rbind(datFinalmol_rep,datt,row.names = NULL)
                    
                        }
                    }else{
                        datt<-cbind(datftemp[qi,1:4],
                        chromosome = chromosome,
                        chromPosition = chromPosition,
                        contigPosition = Position, 
                        ContigID = Contig,
                        contigSiteID = contigSiteID, 
                        datftemp[qi,5:ncol(datftemp)],row.names = NULL)
                        datFinalmol_rep<-rbind(datFinalmol_rep,datt,row.names = NULL)
                
                    }
            ik = ik+1
                }else{
                    #chromosome <- as.numeric(unique(dattemp2$chromosome))
                    datt<-cbind(datftemp[qi,1:4],
                    chromosome = chromo,
                    chromPosition = 0,
                    contigPosition = 0, 
                    ContigID = unique(datmol$ContigID),
                    contigSiteID = 0, 
                    datftemp[qi,5:ncol(datftemp)],row.names = NULL)
                    datFinalmol_rep<-rbind(datFinalmol_rep,datt,row.names = NULL)
            next
            }
            
            }else if(datftemp$LabelChannel[qi] == cntMeth){
                datt<-cbind(datftemp[qi,1:4],
                        chromosome = chromo,
                        chromPosition = 0, 
                        contigPosition = 0,
                        ContigID = unique(datmol$ContigID),                        
                        contigSiteID = 0, 
                        datftemp[qi,5:ncol(datftemp)],row.names = NULL)
                datFinalmol_rep<-rbind(datFinalmol_rep,datt,row.names = NULL)
                
            }else{
                datt<-cbind(datftemp[qi,1:4],
                    chromosome = chromo,
                    chromPosition = 0, 
                    contigPosition = 0,
                    ContigID = unique(datmol$ContigID),                    
                    contigSiteID = 0, 
                    datftemp[qi,5:ncol(datftemp)],row.names = NULL)
                datFinalmol_rep<-rbind(datFinalmol_rep,datt,row.names = NULL)
                
            }
            
        }
        datFinalmol_rep1 <-datFinalmol_rep[order(datFinalmol_rep$SiteID),]
        datFinalmol <- rbind(datFinalmol, datFinalmol_rep1,row.names = NULL)    
        } else if(OrientationMolecule == "+"){
                    for(qi in 1:nrow(datftemp)){
                        #print(qi)
                        if(datftemp$LabelChannel[qi] == cntNick){
                            
                            #dattemp3 <- datftemp[which(as.numeric((datftemp$NickLoc)[qi]) == molsiteID[ik]),]
                            #print(paste("molsiteID[ik]:",molsiteID[ik]))
                            #print(paste("ik",ik))
                            if(is.na(molsiteID[ik])) {
                                molsiteID[ik] = 0
                            } else{
                                molsiteID[ik] = molsiteID[ik]
                            }
                            if(datftemp$NickLoc[qi] == as.character(molsiteID[ik])){
                                dattemp2 <- dattemp1[which((dattemp1$MolsiteID) == molsiteID[ik]),]
                                #print(names(dattemp2))
                                contigSiteID <- as.numeric(dattemp2$ContigSiteID)
                                Position <- (as.numeric(dattemp2$LocRef))
                                Contig <- dattemp2$ContigID
                                #print(Contig)
                                chromPosition <- as.numeric(dattemp2$chrompos)
                                chromosome <- as.numeric(dattemp2$chromosome)
                                if(length(contigSiteID)> 1){
								    'print(1)
								    print(qi)
						            print(OrientationMolecule)'
                                    for(lo in 1:length(contigSiteID)){
                                        datt<-cbind(datftemp[qi,1:4],
                                            chromosome = chromosome[lo],
                                            chromPosition = chromPosition[lo],
                                            contigPosition = Position[lo], 
                                            ContigID = Contig[lo],
                                            contigSiteID = contigSiteID[lo], 
                                            datftemp[qi,5:ncol(datftemp)],row.names = NULL)
                                        datFinalmol<-rbind(datFinalmol,datt,row.names = NULL)
                                        
                                    }
                                }else{
								    'print(2)
								    print(qi)
						            print(OrientationMolecule)'
                                    datt<-cbind(datftemp[qi,1:4],
                                        chromosome = chromosome,
                                        chromPosition = chromPosition,
                                        contigPosition = Position, 
                                        ContigID = Contig,
                                        contigSiteID = contigSiteID, 
                                        datftemp[qi,5:ncol(datftemp)],row.names = NULL)
                                    datFinalmol<-rbind(datFinalmol,datt,row.names = NULL)
                                    
                                }
                                ik = ik+1
                            }else{
							    'print(3)
								print(qi)
						        print(OrientationMolecule)'
                                datt<-cbind(datftemp[qi,1:4],
                                    chromosome = chromo,
                                    chromPosition = 0,
                                    contigPosition = 0, 
                                    ContigID = unique(datmol$ContigID),
                                    contigSiteID = 0, 
                                    datftemp[qi,5:ncol(datftemp)])
                                datFinalmol<-rbind(datFinalmol,datt,row.names = NULL)
                                next
                            }
                                
                        }else if(datftemp$LabelChannel[qi] == cntMeth){
						    'print(4)
							print(qi)
						    print(OrientationMolecule)'
                            datt<-cbind(datftemp[qi,1:4],
                            chromosome = chromo,
                            chromPosition = 0,
                            contigPosition = 0, 
                            ContigID = unique(datmol$ContigID),
                            contigSiteID = 0, 
                            datftemp[qi,5:ncol(datftemp)])
                            datFinalmol<-rbind(datFinalmol,datt,row.names = NULL)
                            
                        }else{
						    'print(5)
							print(qi)
						    print(OrientationMolecule)'
                            datt<-cbind(datftemp[qi,1:4],
                            chromosome = chromo,
                            chromPosition = 0,
                            contigPosition = 0, 
                            ContigID = unique(datmol$ContigID),
                            contigSiteID = 0,
                            datftemp[qi,5:ncol(datftemp)])
                            datFinalmol<-rbind(datFinalmol,datt,row.names = NULL)
                            
                        }
            
    }
                
            }else{
                print("Orientation not present")
            }
            
        }
    #stopImplicitCluster(cl)
	return(datFinalmol)
}
#' Mapping nick label location from contig to reference
#'
#' @param refxmap character. contig reference xmap data.
#' @param refcmap character. reference cmap data.
#' @param contigcmap character. contig cmap data.
#' @param contigID character. query contig ID.
#' @param returnMethod character. Output method. Text or dataframe.
#' @param hmdir integer. Home Directory.
#' @param chrom integer. chromosome number.
#' @return dataframe containing the molecule file, with the reference position
#' @examples
#' refcontigXmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/ContigRef.xmap", package="methometR")
#' refCmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/hg19ref_r.cmap", package="methometR")
#' Contigqcmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/SampContig_q.cmap", package="methometR")
#' refxmapdat <- readingXmap(refcontigXmap)
#' refcmapdat <- readingCmap(refCmap)
#' contigcmapdat <- readingCmap(Contigqcmap)
#' nickRef<-nickReference(refxmap = refxmapdat, refcmap = refcmapdat, 
#'     contigcmap = contigcmapdat, contigID = 6701, 
#'    returnMethod =c("dataFrame"),  chrom = 4)
#' @import utils
#' @import data.table
#' @export
nickReference <- function(refxmap, refcmap, contigcmap, contigID, 
    returnMethod =c("Text", "dataFrame"), hmdir, chrom){
	print("Searching for the nick labels in the reference")
    refxmap = refxmap
    contigcmap <- contigcmap
    contigID = contigID
	chrom = chrom[1]
    #print(contigID)
    #print(contigID)
    dat1 <- refxmap[which(refxmap$QryContigID == contigID
	    & refxmap$RefContigID  == chrom),]
    aln<-dat1$Alignment
    if(length(aln) > 1){
        dat <- data.frame()
        for(z in 1:length(aln)){
            st<-strsplit(as.character(aln[z]),split="\\)\\(")
            #st1<-as.character(st[[1]])
            datt <- data.frame()
            for(i in 1:length(st)){
                st1 <- st[[i]]
                queryloc<-c(); refloc <- c()
                for (ti in 1:length(st1)){
                        if(length(grep("\\(",st1[ti]))==1){
                            st2<-strsplit(st1[ti],split="\\(")
                            st3<-strsplit(st2[[1]][2],split=",")
                            queryloc<-c(queryloc,as.numeric(st3[[1]][2]))
                            refloc<-c(refloc,as.numeric(st3[[1]][1]))
                        }  
                        else if(length(grep("\\)",st1[ti]))==1){
                            st2<-strsplit(st1[ti],split="\\)")
                            st3<-strsplit(st2[[1]][1],split=",")
                            queryloc<-c(queryloc,as.numeric(st3[[1]][2]))
                            refloc<-c(refloc,as.numeric(st3[[1]][1]))
                        }
                        else{
                            st3<-strsplit(st1[ti],split=",")
                            queryloc<-c(queryloc,as.numeric(st3[[1]][2]))
                            refloc<-c(refloc,as.numeric(st3[[1]][1]))
                        }
                    }
                queryLoc = as.numeric(queryloc)
                h<-data.table(queryLoc = as.numeric(queryloc), 
				    refLoc = as.numeric(refloc))
				setkey(h, queryLoc)
                #.set(h, keys=as.numeric(queryloc),values=as.numeric(refloc))
                
                ###Getting unique cmapID 
                contCmapID <- as.numeric(unique(dat1$QryContigID)[i])
                #datFinal <- data.frame()
                datFinalmol <- data.frame()
                LocRef<-c();cmapid<-c();chromosome<-c(); queryref <- c(); siteID <- c()
                k1 <- as.numeric(unique(queryLoc))
                for (ll in 1:length(queryloc)){
                    
                    
                    if(length(grep(paste("^",as.numeric(queryloc[ll]),"$",sep=""),k1))>=1){
                        #print(as.numeric(queryloc[ll]))
                        #val1<-as.numeric(hash::values(h,keys=as.numeric(queryloc[ll])))
                        val1 <- h[.(as.numeric(k1[ll]))]
						ref <- val1$refLoc
                        dattem<-refcmap[which(as.numeric(refcmap$CMapId)==chrom 
                            & as.numeric(refcmap$SiteID) %in% ref),]
                        #print(dim(dattem))
                        LocRef<-c(LocRef,as.numeric(dattem$Position))
                        cmapid<-c(cmapid,as.numeric(contCmapID))
                        chromosome<-c(chromosome,chrom)
                        queryref <- c(queryref, as.numeric(queryloc[ll]))
                    }else{
                    #print(paste("CMapID",CmapID[ui],sep=""))
                        LocRef<-c(LocRef,0)
                        cmapid<-c(cmapid,as.numeric(contCmapID))
                        chromosome<-c(chromosome,chrom)
                        queryref <- c(queryref, as.numeric(queryloc[ll]))
                    }
                }       
                dattemp<-data.frame(Chromosome=chromosome,cmapid=cmapid,LocRef=LocRef, queryref = queryref)
                #write.table(dattemp, paste0("Contig_", unique(cmapid), "_r_mod.cmap"), row.names = FALSE)
                datt <- rbind(datt, dattemp,row.names = NULL)
            }
            dat <- rbind(dat,datt,row.names = NULL)
        }
    }else{
        st<-strsplit(as.character(aln),split="\\)\\(")
            #st1<-as.character(st[[1]])
            dat <- data.frame()
            for(i in 1:length(st)){
                st1 <- st[[i]]
                queryloc<-c(); refloc <- c()
                for (ti in 1:length(st1)){
                        if(length(grep("\\(",st1[ti]))==1){
                            st2<-strsplit(st1[ti],split="\\(")
                            st3<-strsplit(st2[[1]][2],split=",")
                            queryloc<-c(queryloc,as.numeric(st3[[1]][2]))
                            refloc<-c(refloc,as.numeric(st3[[1]][1]))
                        }  
                        else if(length(grep("\\)",st1[ti]))==1){
                            st2<-strsplit(st1[ti],split="\\)")
                            st3<-strsplit(st2[[1]][1],split=",")
                            queryloc<-c(queryloc,as.numeric(st3[[1]][2]))
                            refloc<-c(refloc,as.numeric(st3[[1]][1]))
                        }
                        else{
                            st3<-strsplit(st1[ti],split=",")
                            queryloc<-c(queryloc,as.numeric(st3[[1]][2]))
                            refloc<-c(refloc,as.numeric(st3[[1]][1]))
                        }
                    }
            
                'h<-hash()
                .set(h, keys=as.numeric(queryloc),values=as.numeric(refloc))'
				queryLoc = as.numeric(queryloc)
				h<-data.table(queryLoc = as.numeric(queryloc), 
				    refLoc = as.numeric(refloc))
				setkey(h, queryLoc)
                
                ###Getting unique cmapID 
                contCmapID <- as.numeric(unique(dat1$QryContigID)[i])
                #datFinal <- data.frame()
                datFinalmol <- data.frame()
                LocRef<-c();cmapid<-c();chromosome<-c(); queryref <- c(); siteID <- c()
                #k1<-as.numeric(hash::keys(h))
				k1 <- as.numeric(unique(h$queryLoc))
                for (ll in 1:length(queryloc)){
                    
                    
                    if(length(grep(paste("^",as.numeric(queryloc[ll]),"$",sep=""),k1))>=1){
                        #print(as.numeric(queryloc[ll]))
                        val1 <- h[.(as.numeric(k1[ll]))]
						ref <- val1$refLoc
                        dattem<-refcmap[which(as.numeric(refcmap$CMapId)==chrom 
                            & as.numeric(refcmap$SiteID) %in% ref),]
                        #print(dim(dattem))
                        LocRef<-c(LocRef,as.numeric(dattem$Position))
                        cmapid<-c(cmapid,as.numeric(contCmapID))
                        chromosome<-c(chromosome,chrom)
                        queryref <- c(queryref, as.numeric(queryloc[ll]))
                    }else{
                    #print(paste("CMapID",CmapID[ui],sep=""))
                        LocRef<-c(LocRef,0)
                        cmapid<-c(cmapid,as.numeric(contCmapID))
                        chromosome<-c(chromosome,chrom)
                        queryref <- c(queryref, as.numeric(queryloc[ll]))
                    }
                }
                dattemp<-data.frame(Chromosome=chromosome,cmapid=cmapid,LocRef=LocRef, queryref = queryref)
                #write.table(dattemp, paste0("Contig_", unique(cmapid), "_r_mod.cmap"), row.names = FALSE)
                dat <- rbind(dat, dattemp,row.names = NULL)
            }
    }
    ####Methylation generation
    cmapID<-as.numeric(unique(dat$cmapid))
    datFinal <- data.frame()
    #dattemp1<-dat[which(as.numeric(dat$cmapid)==cmapID),]
    
    if(length(cmapID) > 1){
        for(oo in 1:length(cmapID)){
        #print(oo)
        datftemp <- contigcmap[which(as.numeric(contigcmap$CMapId) == cmapID[oo]),]
        #chromosome=dattemp1$Chromosome
        datc <- dat[which(dat$cmapid ==cmapID[oo]),]
        Position=(as.numeric(datc$LocRef))
        qry <- (as.numeric(datc$queryref))
        pqry <- paste0("^", qry , "$")
            for( yy in 1:nrow(datftemp)){
                if(as.numeric(datftemp$LabelChannel[yy])==1){
                        passite <- paste0("^", datftemp$SiteID[yy], "$")
                        g1 <- grep(passite, pqry, fixed = TRUE)
                            if(length(g1) == 1){
                                datt<-cbind(datftemp[yy,1:6],
                                    chromosome = chrom,
                                    chromosomePosition = Position[g1],
                                    datftemp[yy,7:ncol(datftemp)])
                                datFinal<-rbind(datFinal,datt, row.names = NULL)
                            }
                            else if(length(g1) > 1){
                                dattt <- data.frame()
                                for (nn in 1:length(g1)){
                                    datt<-cbind(datftemp[yy,1:6],
                                        chromosome = chrom,
                                        chromosomePosition = Position[g1[nn]],
                                        datftemp[yy,7:ncol(datftemp)])
                                    dattt <- rbind(dattt, datt, row.names = NULL)
                                }
                                 datFinal<-rbind(datFinal,dattt, row.names = NULL)
                            }else {
                              
                              'print(paste0("Locations dont match :", 
                                datftemp$SiteID[yy]))'
                                datt<-cbind(datftemp[yy,1:6],
                                        chromosome = chrom,
                                        chromosomePosition = NA,
                                        datftemp[yy,7:ncol(datftemp)])
                                datFinal<-rbind(datFinal,datt, row.names = NULL)
                          }
                    }
                else if(as.numeric(datftemp$LabelChannel[yy]) == 0){
                    posn <- Position[yy-1]
                    pos <- posn + (datftemp$Position[yy]-datftemp$Position[yy-1])
                    datt<-cbind(datftemp[yy,1:6],
                        chromosome=chrom,
                        chromosomePosition=pos,
                        datftemp[yy,7:ncol(datftemp)])
                    datFinal<-rbind(datFinal,datt, row.names = NULL)
                }
                else{
                    print("Inappropriate label")
                }
            }
            'write.table(datFinal, file.path(hmdir, paste0("Contigs_",
                cmapID[oo], "_q.cmap")), 
                sep = "\t",row.names = FALSE)'
        }
    
    } else if (length(cmapID) == 1){
        datftemp <- contigcmap[which(as.numeric(contigcmap$CMapId)==cmapID),]
        #chromosome=dattemp1$Chromosome
        datc <- dat[which(dat$cmapid ==cmapID),]
        Position=(as.numeric(datc$LocRef))
        qry <- (as.numeric(datc$queryref))
        pqry <- paste0("^", qry , "$")
            for( yy in 1:nrow(datftemp)){
			    'print(cmapID)
			    print(yy)'
                if(as.numeric(datftemp$LabelChannel[yy])==1){
                        passite <- paste0("^", datftemp$SiteID[yy], "$")
                        g1 <- grep(passite, pqry, fixed = TRUE)
                            if(length(g1) == 1){
                                datt<-cbind(datftemp[yy,1:6],
                                  chromosome = chrom,
                                  chromosomePosition = Position[g1],
                                  datftemp[yy,7:ncol(datftemp)])
                                datFinal<-rbind(datFinal,datt, row.names = NULL)
                            }
                            else if(length(g1) > 1){
						    for(pp in 1:length(g1)){
                                datt<-cbind(datftemp[yy,1:6],
                                    chromosome = chrom,
                                    chromosomePosition = Position[g1[pp]],
                                    datftemp[yy,7:ncol(datftemp)])
                            datFinal<-rbind(datFinal,datt, row.names = NULL)
							}
                            }else {
                                
                                #print(paste0("Locations don't match :", 
                               #  datftemp$SiteID[yy]))
                                datt<-cbind(datftemp[yy,1:6],
                                      chromosome = chrom,
                                      chromosomePosition = NA,
                                      datftemp[yy,7:ncol(datftemp)])
                                 datFinal<-rbind(datFinal,datt, row.names = NULL)
                            }
                    }
                else if(as.numeric(datftemp$LabelChannel[yy]) == 0){
                    posn <- Position[yy-1]
                    pos <- posn + (datftemp$Position[yy]-datftemp$Position[yy-1])
                    datt<-cbind(datftemp[yy,1:6],
                        chromosome=chrom,
                        chromosomePosition=pos,
                        datftemp[yy,7:ncol(datftemp)])
                    datFinal<-rbind(datFinal,datt, row.names = NULL)
                }
                else{
                    print("Inappropriate label")
                }
            }
       
    } else{ print("No cmap ID")}
    if(returnMethod == "Text"){
        write.table(datFinal, file.path(hmdir, paste0("Contigs_",
            cmapID, "_q.cmap")), 
            sep = "\t",row.names = FALSE)
    }else if(returnMethod == "dataFrame"){
        return(datFinal)
    }else{stop("No returnMethodSpecified")}
#return(datFinal)
}
#' Liftover of molecule methylation and nick label position to contig coordinates
#'
#' @param molcmap character. molecule cmap data.
#' @param outputType character. Output method. dataframe or text.
#' @param cntNick integer. Nick label value.
#' @param cntMeth integer. Meth label value.
#' @return dataframe containing the molecule file, with the reference position
#' @examples
#' xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
#' cmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
#' modcmap <- readingCmap(cmap)
#' modxmap <- readingXmap(xmap)
#' modMolcmap <- modmolcmap(molcmap = modcmap,
#'   xmapdata = modxmap,
#'   cntNick = 1,
#'   cntMeth=2)
#' refcontigXmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/ContigRef.xmap", package="methometR")
#' refCmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/hg19ref_r.cmap", package="methometR")
#' Contigqcmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/SampContig_q.cmap", package="methometR")
#' refxmapdat <- readingXmap(refcontigXmap)
#' refcmapdat <- readingCmap(refCmap)
#' contigcmapdat <- readingCmap(Contigqcmap)
#' nickRef<-nickReference(refxmap = refxmapdat, refcmap = refcmapdat, 
#'     contigcmap = contigcmapdat, contigID = 6701, 
#'    returnMethod =c("dataFrame"),  chrom = 4)
#' xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
#' molcmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
#' contigcmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContig_r.cmap", package="methometR")
#' contigID = 6701; cntNick = 1; cntMeth=2
#' contigmodcmap <- readingXmap(contigcmap)
#' modmolcmapNick <-molNickLoc(nickrefLoc = nickRef,
#'    molcmap = modMolcmap,
#'    xmapdata = modxmap, 
#'    contigID = contigID, cntNick = cntNick, 
#'    cntMeth = cntMeth)
#' MethylContigLoc<-molMethylLoc(molcmap = modmolcmapNick,
#'    outputType = c("dataframe"),
#'    cntNick = cntNick, 
#'    cntMeth = cntMeth)
#' @import utils
#' @export
molMethylLoc <- function(molcmap, outputType = c("dataframe", "text"), cntMeth, cntNick){
    print("Searching for corresponding methylation label locations in the contig")
    molcmap = molcmap
    dataContig <- data.frame()
    datFinalmol = molcmap
    cmapID<-as.character(unique(datFinalmol$CMapId))
    #print(paste("cmapID::",cmapID))
    ###Adding the probable location of Methylation
    #datcontig<-data.frame()
	'cores=detectCores()
    cl <- makeCluster(5) #not to overload your computer
	print(cl)
    registerDoParallel(cl)'
	#dataContig <- foreach(ti=1:length(cmapID), .combine=rbind) %dopar% {
        for ( ti in 1:length(cmapID)){
            'print(paste("ti", ti))
	    	print(paste("cmapID", cmapID[ti]))'
            datcontig<-data.frame()
            datcontig<-datFinalmol[which(as.numeric(datFinalmol$CMapId) == as.numeric(cmapID[ti])),]
            orientMolecule <- unique(as.character(datcontig$OrientationMolecule))
            #print(orientMolecule)
            datcontig[is.na(datcontig)]<-0
            
            for( y in 1:nrow(datcontig)){
                z=y-1
                x=y+1
                #print(paste("y", y))
                if(y == 1){
	    		    
                    if((datcontig$LabelChannel[y] == cntMeth 
	    			    &  datcontig$contigPosition[y] == 0)
                        & (datcontig$LabelChannel[x] == cntNick 
                        & datcontig$contigPosition[x]>0)){
                        if(orientMolecule == "+"){
                            datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                            if(datcontig$contigPosition[y] > datcontig$contigPosition[x]| datcontig$contigPosition[y] < 0){
                                datcontig$contigPosition[y] = 0
                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                        }else{
                            datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                            if(datcontig$contigPosition[y] < 0){
                                datcontig$contigPosition[y] = 0
                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                        
                        }
                    }
                    else if((datcontig$LabelChannel[y]== cntMeth 
                        & datcontig$contigPosition[y] == 0) 
                        & (datcontig$LabelChannel[x] == cntNick 
                        & datcontig$contigPosition[x] == 0)){
                        g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(y < x & datcontig$contigPosition[x] > 0){
                                if(orientMolecule == "+"){
                                    datcontig$contigPosition[y]=as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$contigPosition[y] > datcontig$contigPosition[x]| datcontig$contigPosition[y] < 0){
                                        datcontig$contigPosition[y] = 0
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                            } else{
                                    datcontig$contigPosition[y]=as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$contigPosition[y] < 0){
                                        datcontig$contigPosition[y] = 0
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                }
                                
                                ji=length(g1)+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }            
                    
                    }
                    else if(datcontig$LabelChannel[y]== cntMeth & (datcontig$LabelChannel[x] == cntMeth & datcontig$contigPosition[x]==0)){
                        g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            ##print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(y < x & datcontig$contigPosition[x] > 0){
                                if((orientMolecule == "+")){
                                    datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                    if(datcontig$contigPosition[y] > datcontig$contigPosition[x]| datcontig$contigPosition[y] < 0){
                                        datcontig$contigPosition[y] = 0
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                } else{
                                    datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                    if(datcontig$contigPosition[y] < 0){
                                        datcontig$contigPosition[y] = 0
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                                }
                            ji=length(g1)+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }
                    }
                    if(datcontig$LabelChannel[y] == cntNick 
                        & (datcontig$LabelChannel[x] == cntNick 
                        & datcontig$contigPosition[x]>0)){
                        if(orientMolecule == "+"){
                            datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                            if(datcontig$contigPosition[y] > datcontig$contigPosition[x]| datcontig$contigPosition[y] < 0){
                                datcontig$contigPosition[y] = 0
                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                        }else{
                            datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                            if(datcontig$contigPosition[y] < 0){
                                datcontig$contigPosition[y] = 0
                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                        
                        }
                    }
                    else if((datcontig$LabelChannel[y]== cntNick 
                        & datcontig$contigPosition[y] == 0) 
                        & (datcontig$LabelChannel[x] == cntNick 
                        & datcontig$contigPosition[x] == 0)){
                        g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(y < x & datcontig$contigPosition[x] > 0){
                                if(orientMolecule == "+"){
                                    datcontig$contigPosition[y]=as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$contigPosition[y] > datcontig$contigPosition[x]| datcontig$contigPosition[y] < 0){
                                        datcontig$contigPosition[y] = 0
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                            } else{
                                    datcontig$contigPosition[y]=as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$contigPosition[y] < 0){
                                        datcontig$contigPosition[y] = 0
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                }
                                
                                ji=length(g1)+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }            
                    
                    }
                    else if(datcontig$LabelChannel[y]== cntNick 
                        & (datcontig$LabelChannel[x] == cntMeth 
                        & datcontig$contigPosition[x]==0)){
                        g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            ##print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(y < x & datcontig$contigPosition[x] > 0){
                                if((orientMolecule == "+")){
                                    datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                    if(datcontig$contigPosition[y] > datcontig$contigPosition[x]| datcontig$contigPosition[y] < 0){
                                        datcontig$contigPosition[y] = 0
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                } else{
                                    datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                    if(datcontig$contigPosition[y] < 0){
                                        datcontig$contigPosition[y] = 0
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                                }
                            ji=length(g1)+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }
                    }
                    else if (datcontig$LabelChannel[y] == cntNick & datcontig$contigPosition[y]>0){
                        datcontig$contigPosition[y] = datcontig$contigPosition[y]
                    }
                    else {
                    
                        datcontig$contigPosition[y]=0
                    }
                }
                else{
                    if(datcontig$LabelChannel[y] == cntMeth 
                        & (datcontig$LabelChannel[z] == cntNick 
                        & datcontig$contigPosition[z]>0)){
                    
                        'if(orientMolecule == "+"){'
                            'datcontig$contigLeft[y] <- as.numeric(datcontig$contigPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))'
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:nrow(datcontig)]))
                            posn <- as.numeric(datcontig$Position[y+1:nrow(datcontig)])
                            contposn <- as.numeric(datcontig$contigPosition[y+1:nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[x] > 0  
                                    & (datcontig$contigPosition[y]== 0 
                                    | is.nan(datcontig$contigPosition[y]))){
                                    datcontig$contigPosition[y] <- (
                        (
                (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))/(as.numeric(posn[x]) - as.numeric(datcontig$Position[z])))* (contposn[x] - datcontig$contigPosition[z]))+ datcontig$contigPosition[z]
                                ji=length(g1)+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }
                            if(datcontig$contigPosition[y]== 0 
                                | is.nan(datcontig$contigPosition[y])){
                                if((orientMolecule == "+")){
                                       datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                       if(datcontig$contigPosition[y] < datcontig$contigPosition[z]| datcontig$contigPosition[y] < 0){
                                           datcontig$contigPosition[y] = 0
                                       }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                  } else{
                                       datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$contigPosition[y] < 0){
                                              datcontig$contigPosition[y] = 0
                                          }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                                   }
                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                            'if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                datcontig$contigPosition[y] <- datcontig$contigLeft[y]
                            }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                datcontig$contigPosition[y] <- datcontig$contigRight[y]
                            }
                            else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                datcontig$contigRight[y] = datcontig$contigRight[y]
                            }
                            else{
                                datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                            }'
                        
                               '} else{
                                datcontig$contigLeft[y] = as.numeric(datcontig$contigPosition[z]) - (as.numeric(datcontig$Position[z])-as.numeric(datcontig$Position[y]))
                                g1<-grep(cntNick, as.character(datcontig$LabelChannel[(y+1):length(datcontig$LabelChannel)]))
                                posn <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                ji=1
                                while(ji <=length(g1)){
                                    x=g1[ji]
                                    #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(contposn[x]>0  & datcontig$contigPosition[y] == 0){
                                        if(orientMolecule == "+"){
                                            datcontig$contigRight[y] = as.numeric(contposn[x]) - (as.numeric(posn[x])-as.numeric(datcontig$Position[y]))
                                            datcontig$contigPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                         } else{
                                             datcontig$contigRight[y] = as.numeric(contposn[x]) + (as.numeric(datcontig$Position[y])-as.numeric(posn[x]))
                                            datcontig$contigPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                        }
                                    ji=length(g1)+1
                                    }
                                    else{
                                        ji=ji+1
                                        next
                                    } 
                                }
                                
                                if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                    datcontig$contigPosition[y] <- datcontig$contigLeft[y]
                                }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigPosition[y] <- datcontig$contigRight[y]
                                }
                                else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                    datcontig$contigRight[y] = datcontig$contigRight[y]
                                }
                                 else{
                                    datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                                }
                            }'
                    }
                    else if(datcontig$LabelChannel[y] == cntMeth & 
                        (datcontig$LabelChannel[x] == cntNick 
                        & datcontig$contigPosition[x]>0)){
                        
                            g1<-grep(cntNick, as.character(datcontig$LabelChannel[1:z]))
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                                ji= length(g1)
                                while(ji > 0){
                                    z=g1[ji]
                                    #print (paste(z,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(contposn[z] > 0 
                                        & (datcontig$contigPosition[y]== 0 
                                        | is.nan(datcontig$contigPosition[y]))){
                                         datcontig$contigPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(datcontig$Position[x]) - as.numeric(posn[z])))* (datcontig$contigPosition[x] - contposn[z])) + contposn[z]
                                        ji= -1
                                    }
                                    else{
                                        ji=ji-1
                                        next
                                    }
                                }
                                if(datcontig$contigPosition[y]== 0 
                                | is.nan(datcontig$contigPosition[y])){
                                    if((orientMolecule == "+")){
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                           if(datcontig$contigPosition[y] > datcontig$contigPosition[x]
                                               | datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                      } else{
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                               if(datcontig$contigPosition[y] < 0){
                                                  datcontig$contigPosition[y] = 0
                                              }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    
                                       }
                                }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                            'if(orientMolecule == "+"){
                                #datcontig$contigRight[y]=as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                
                                
                                
                                if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                    datcontig$contigPosition[y] <- datcontig$contigLeft[y]
                                }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigPosition[y] <- datcontig$contigRight[y]
                                }
                                else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                    datcontig$contigRight[y] = datcontig$contigRight[y]
                                }
                                 else{
                                    datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                                }
                                } else{
                                    datcontig$contigRight[y] = as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[y])-as.numeric(datcontig$Position[x]))
                                    g1<-grep("1",as.character(datcontig$LabelChannel[1:z]))
                                    ji= length(g1)
                                    posn <- as.numeric(datcontig$Position[1:z])
                                    contposn <- as.numeric(datcontig$contigPosition[1:z])
                                    while(ji > 0){
                                        z=g1[ji]
                                        #print (paste(z,":",y,":",cmapID[ti],":",ji,sep=""))
                                        if(contposn[z]>0  & datcontig$contigPosition[y] == 0){
                                            if(orientMolecule == "+"){
                                                datcontig$contigLeft[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y])-as.numeric(posn[z]))
                                                datcontig$contigPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                             } else{
                                                 datcontig$contigLeft[y] = as.numeric(contposn[z]) - (as.numeric(posn[z])-as.numeric(datcontig$Position[y]))
                                                datcontig$contigPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                            }
                                            ji = -1
                                        }
                                        else{
                                            ji=ji - 1
                                            next
                                        }
                                    }
                                    if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                        datcontig$contigPosition[y] <- datcontig$contigLeft[y]
                                    }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                        datcontig$contigPosition[y] <- datcontig$contigRight[y]
                                    }
                                    else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                        datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                        datcontig$contigRight[y] = datcontig$contigRight[y]
                                    }
                                     else{
                                        datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                                    }
                                }       
                    '
                    }
                    else if(datcontig$LabelChannel[y] == cntMeth 
                        & ((datcontig$LabelChannel[z] == cntNick 
                        | datcontig$LabelChannel[x] == cntNick) 
                        & (datcontig$contigPosition[z] == 0
                        |datcontig$contigPosition[x]==0))){
                        g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                        g2 <- grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                        if((datcontig$contigPosition[y]== 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g1) >=1 & length(g2) >=1)){
                            ji=length(g1)
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                            while(ji > 0){
                                z=g1[ji]
                                
                                #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[z] > 0 
                                    & (datcontig$contigPosition[y]== 0 
                                    | is.nan(datcontig$contigPosition[y]))){
                                        #print(ji)
                                        #datcontig$contigLeft[y]=as.numeric(contposn[z]) + (as.numeric(posn[z]) - as.numeric(datcontig$Position[y]))
                                        ###Calculating x
                                        g11<-grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                                        posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                        contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                        ki=1
                                        while(ki <=length(g11)){
                                            x1=g11[ki]
                                            ##print (paste(x1,":",y,":",cmapID[ti],":",ki,sep=""))
                                            if(contposn1[x1] > 0 
                                                & (datcontig$contigPosition[y]== 0 
                                                | is.nan(datcontig$contigPosition[y]))){
                                                datcontig$contigPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(posn1[x1]) - as.numeric(posn[z])))* (contposn1[x1] - contposn[z]))+ contposn[z]
                                                ki=length(g11)+1
                                                ji =  -1
                                            }
                                            else{
                                                ki=ki+1
                                                next
                                            }
                                        }
                                    if(datcontig$contigPosition[y]== 0 
                                        | is.nan(datcontig$contigPosition[y])){
                                        if((orientMolecule == "+")){
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$contigPosition[y] < datcontig$contigPosition[z]| datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        } else{
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                               if(datcontig$contigPosition[y] < 0){
                                                  datcontig$contigPosition[y] = 0
                                              }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        
                                        }
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}    
                                     
                                    
                                    ji =   ji - 1
                                    #print(paste("1:",ji))
                                }
                                else{
                                    datcontig$contigPosition[y]=0
                                    ji = ji-1
                                    #print(paste("2:",ji))
                                }
                            }
                            if(datcontig$contigPosition[y]== 0 
                                | is.nan(datcontig$contigPosition[y])){
                                ji=length(g1)
                                g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                                posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                ji=1
                                while(ji <=length(g1)){
                                    x=g1[ji]
                                    #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(y < x & datcontig$contigPosition[x] > 0){
                                        if(orientMolecule == "+"){
                                            datcontig$contigPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                            if(datcontig$contigPosition[y] > contposn1[x] 
                                                | datcontig$contigPosition[y] < 0){
                                                datcontig$contigPosition[y] = 0
                                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    } else{
                                            datcontig$contigPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                            if(datcontig$contigPosition[y] < 0){
                                                datcontig$contigPosition[y] = 0
                                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        }
                                        
                                        ji = ji+1
                                    }
                                    else{
                                        ji=ji+1
                                        next
                                    }
                                }
                            }else{
                                datcontig$contigPosition[y] = datcontig$contigPosition[y]
                            }    
                            
                        } else if((datcontig$contigPosition[y] == 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g2) >=1 & length(g1) == 0)){
                                ji=length(g1)
                                g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                                posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y < x & datcontig$contigPosition[x] > 0){
                                    if(orientMolecule == "+"){
                                        datcontig$contigPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$contigPosition[y] > contposn1[x] 
                                            | datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                } else{
                                        datcontig$contigPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    }
                                    
                                    ji = ji+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }            
            
                        } else if((datcontig$contigPosition[y]== 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g2) == 0 & length(g1) >= 1)){
                            g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                            ji=length(g1)
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                            while(ji > 0){
                                z=g1[ji]
                                #print(ji)
                                #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y > z & (contposn[z] > 0 & !is.na(contposn[z]))){
                                    if((orientMolecule == "+")){
                                        datcontig$contigPosition[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                        if(datcontig$contigPosition[y] < contposn[z]| datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                   } else{
                                        datcontig$contigPosition[y] = as.numeric(contposn[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                            if(datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                                    }
                                       ji= ji - 1
                                } else{
                                    ji = ji - 1
                                    next
                                }               
                            }
                        }else{
                            #print("Locations do not match!!!!!")
                        }    
                    }
                    else if(datcontig$LabelChannel[y] == cntMeth & 
                        ((datcontig$LabelChannel[z] == cntMeth 
                        | datcontig$LabelChannel[x]== cntMeth) 
                        & (datcontig$contigPosition[z]==0
                        |datcontig$contigPosition[x]==0))){
                        g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                        g2 <- grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                        if((datcontig$contigPosition[y]== 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g1) >=1 & length(g2) >=1)){
                            ji=length(g1)
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                            while(ji > 0){
                                z=g1[ji]
                                #print(ji)
                                #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[z] > 0 
                                    & (datcontig$contigPosition[y]== 0 
                                    | is.nan(datcontig$contigPosition[y]))){
                                    
                                        #datcontig$contigLeft[y]=as.numeric(contposn[z]) + (as.numeric(posn[z]) - as.numeric(datcontig$Position[y]))
                                        ###Calculating x
                                        g11<-grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                                        posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                        contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                        ki=1
                                        while(ki <=length(g11)){
                                            x1=g11[ki]
                                            ##print (paste(x1,":",y,":",cmapID[ti],":",ki,sep=""))
                                            if(contposn1[x1] > 0 
                                                & (datcontig$contigPosition[y]== 0 
                                                | is.nan(datcontig$contigPosition[y]))){
                                                datcontig$contigPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(posn1[x1]) - as.numeric(posn[z])))* (contposn1[x1] - contposn[z]))+ contposn[z]
                                                ki=length(g11)+1
                                                ji = -1
                                            }
                                            else{
                                                ki=ki+1
                                                next
                                            }
                                        }
                                    if(datcontig$contigPosition[y]== 0 
                                        | is.nan(datcontig$contigPosition[y])){
                                        if((orientMolecule == "+")){
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$contigPosition[y] < datcontig$contigPosition[z]| datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        } else{
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                               if(datcontig$contigPosition[y] < 0){
                                                  datcontig$contigPosition[y] = 0
                                              }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        
                                        }
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}    
                                     
                                    
                                    ji =  ji - 1
                                    #print(paste("1:",ji))
                                }
                                else{
                                    datcontig$contigPosition[y]=0
                                    ji = ji-1
                                   # print(paste("2:",ji))
                                }
                            }
                            if(datcontig$contigPosition[y]== 0 
                                | is.nan(datcontig$contigPosition[y])){
                                ji=length(g1)
                                g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                                posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                ji=1
                                while(ji <=length(g1)){
                                    x=g1[ji]
                                    #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(y < x & datcontig$contigPosition[x] > 0){
                                        if(orientMolecule == "+"){
                                            datcontig$contigPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                            if(datcontig$contigPosition[y] > contposn1[x] 
                                                | datcontig$contigPosition[y] < 0){
                                                datcontig$contigPosition[y] = 0
                                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    } else{
                                            datcontig$contigPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                            if(datcontig$contigPosition[y] < 0){
                                                datcontig$contigPosition[y] = 0
                                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        }
                                        
                                        ji = ji+1
                                    }
                                    else{
                                        ji=ji+1
                                        next
                                    }
                                }
                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}    
                            
                        } else if((datcontig$contigPosition[y] == 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g2) >=1 & length(g1) == 0)){
                                ji=length(g1)
                                g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                                posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y < x & datcontig$contigPosition[x] > 0){
                                    if(orientMolecule == "+"){
                                        datcontig$contigPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$contigPosition[y] > contposn1[x] 
                                            | datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                } else{
                                        datcontig$contigPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    }
                                    
                                    ji = ji+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }            
            
                        } else if((datcontig$contigPosition[y]== 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g2) == 0 & length(g1) >= 1)){
                            g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                            ji=length(g1)
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                            while(ji > 0){
                                z=g1[ji]
                                #print(ji)
                                #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y > z & (contposn[z] > 0 & !is.na(contposn[z]))){
                                    if((orientMolecule == "+")){
                                        datcontig$contigPosition[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                        if(datcontig$contigPosition[y] < contposn[z]| datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                   } else{
                                        datcontig$contigPosition[y] = as.numeric(contposn[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                            if(datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                                    }
                                       ji= ji - 1
                                } else{
                                    ji = ji - 1
                                    next
                                }               
                            }
                        }else{
                            #print("Locations do not match!!!!!")
                        }    
                    } else if(datcontig$LabelChannel[y] == cntNick
                        & ((datcontig$LabelChannel[z] == cntNick 
                        | datcontig$LabelChannel[x] == cntNick) 
                        & (datcontig$contigPosition[z] == 0
                        |datcontig$contigPosition[x]==0))){
                        g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                        g2 <- grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                        if((datcontig$contigPosition[y]== 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g1) >=1 & length(g2) >=1)){
                            ji=length(g1)
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                            while(ji > 0){
                                z=g1[ji]
                                
                                #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[z] > 0 
                                    & (datcontig$contigPosition[y]== 0 
                                    | is.nan(datcontig$contigPosition[y]))){
                                    
                                        #datcontig$contigLeft[y]=as.numeric(contposn[z]) + (as.numeric(posn[z]) - as.numeric(datcontig$Position[y]))
                                        ###Calculating x
                                        g11<-grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                                        posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                        contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                        ki=1
                                        while(ki <=length(g11)){
                                            x1=g11[ki]
                                            ##print (paste(x1,":",y,":",cmapID[ti],":",ki,sep=""))
                                            if(contposn1[x1] > 0 
                                                & (datcontig$contigPosition[y]== 0 
                                                | is.nan(datcontig$contigPosition[y]))){
                                                datcontig$contigPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(posn1[x1]) - as.numeric(posn[z])))* (contposn1[x1] - contposn[z]))+ contposn[z]
                                                ki=length(g11)+1
                                                ji =  -1
                                            }
                                            else{
                                                ki=ki+1
                                                next
                                            }
                                        }
                                    if(datcontig$contigPosition[y]== 0 
                                        | is.nan(datcontig$contigPosition[y])){
                                        if((orientMolecule == "+")){
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$contigPosition[y] < datcontig$contigPosition[z]| datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        } else{
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                               if(datcontig$contigPosition[y] < 0){
                                                  datcontig$contigPosition[y] = 0
                                              }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        
                                        }
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}    
                                     
                                    
                                    ji =  ji - 1
                                    #print(paste("1:",ji))
                                }
                                else{
                                    datcontig$contigPosition[y]=0
                                    ji = ji-1
                                    #print(paste("2:",ji))
                                }
                            }
                            if(datcontig$contigPosition[y]== 0 
                                | is.nan(datcontig$contigPosition[y])){
                                ji=length(g1)
                                g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                                posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                ji=1
                                while(ji <=length(g1)){
                                    x=g1[ji]
                                    #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(y < x & datcontig$contigPosition[x] > 0){
                                        if(orientMolecule == "+"){
                                            datcontig$contigPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                            if(datcontig$contigPosition[y] > contposn1[x] 
                                                | datcontig$contigPosition[y] < 0){
                                                datcontig$contigPosition[y] = 0
                                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    } else{
                                            datcontig$contigPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                            if(datcontig$contigPosition[y] < 0){
                                                datcontig$contigPosition[y] = 0
                                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        }
                                        
                                        ji = ji+1
                                    }
                                    else{
                                        ji=ji+1
                                        next
                                    }
                                }
                            }else{
                                datcontig$contigPosition[y] = datcontig$contigPosition[y]
                            }    
                            
                        } else if((datcontig$contigPosition[y] == 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g2) >=1 & length(g1) == 0)){
                                ji=length(g1)
                                g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                                posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y < x & datcontig$contigPosition[x] > 0){
                                    if(orientMolecule == "+"){
                                        datcontig$contigPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$contigPosition[y] > contposn1[x] 
                                            | datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                } else{
                                        datcontig$contigPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    }
                                    
                                    ji = ji+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }            
            
                        } else if((datcontig$contigPosition[y]== 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g2) == 0 & length(g1) >= 1)){
                            g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                            ji=length(g1)
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                            while(ji > 0){
                                z=g1[ji]
                                #print(ji)
                                #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y > z & (contposn[z] > 0 & !is.na(contposn[z]))){
                                    if((orientMolecule == "+")){
                                        datcontig$contigPosition[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                        if(datcontig$contigPosition[y] < contposn[z]| datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                   } else{
                                        datcontig$contigPosition[y] = as.numeric(contposn[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                            if(datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                                    }
                                       ji= ji - 1
                                } else{
                                    ji = ji - 1
                                    next
                                }               
                            }
                        }else{
                            #print("Locations do not match!!!!!")
                        }    
                    }
                    else if(datcontig$LabelChannel[y] == cntNick 
                        & (datcontig$LabelChannel[z] == cntNick 
                        & datcontig$contigPosition[z]>0)){
                    
                        'if(orientMolecule == "+"){'
                            'datcontig$contigLeft[y] <- as.numeric(datcontig$contigPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))'
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:nrow(datcontig)]))
                            posn <- as.numeric(datcontig$Position[y+1:nrow(datcontig)])
                            contposn <- as.numeric(datcontig$contigPosition[y+1:nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[x] > 0  
                                    & (datcontig$contigPosition[y]== 0 
                                    | is.nan(datcontig$contigPosition[y]))){
                                    datcontig$contigPosition[y] <- (
                        (
                (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))/(as.numeric(posn[x]) - as.numeric(datcontig$Position[z])))* (contposn[x] - datcontig$contigPosition[z]))+ datcontig$contigPosition[z]
                                ji=length(g1)+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }
                            if(datcontig$contigPosition[y]== 0 
                                | is.nan(datcontig$contigPosition[y])){
                                if((orientMolecule == "+")){
                                       datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                       if(datcontig$contigPosition[y] < datcontig$contigPosition[z]| datcontig$contigPosition[y] < 0){
                                           datcontig$contigPosition[y] = 0
                                       }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                  } else{
                                       datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$contigPosition[y] < 0){
                                              datcontig$contigPosition[y] = 0
                                          }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                                   }
                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                            'if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                datcontig$contigPosition[y] <- datcontig$contigLeft[y]
                            }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                datcontig$contigPosition[y] <- datcontig$contigRight[y]
                            }
                            else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                datcontig$contigRight[y] = datcontig$contigRight[y]
                            }
                            else{
                                datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                            }'
                        
                               '} else{
                                datcontig$contigLeft[y] = as.numeric(datcontig$contigPosition[z]) - (as.numeric(datcontig$Position[z])-as.numeric(datcontig$Position[y]))
                                g1<-grep(cntNick, as.character(datcontig$LabelChannel[(y+1):length(datcontig$LabelChannel)]))
                                posn <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                ji=1
                                while(ji <=length(g1)){
                                    x=g1[ji]
                                    #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(contposn[x]>0  & datcontig$contigPosition[y] == 0){
                                        if(orientMolecule == "+"){
                                            datcontig$contigRight[y] = as.numeric(contposn[x]) - (as.numeric(posn[x])-as.numeric(datcontig$Position[y]))
                                            datcontig$contigPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                         } else{
                                             datcontig$contigRight[y] = as.numeric(contposn[x]) + (as.numeric(datcontig$Position[y])-as.numeric(posn[x]))
                                            datcontig$contigPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                        }
                                    ji=length(g1)+1
                                    }
                                    else{
                                        ji=ji+1
                                        next
                                    } 
                                }
                                
                                if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                    datcontig$contigPosition[y] <- datcontig$contigLeft[y]
                                }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigPosition[y] <- datcontig$contigRight[y]
                                }
                                else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                    datcontig$contigRight[y] = datcontig$contigRight[y]
                                }
                                 else{
                                    datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                                }
                            }'
                    }
                    else if(datcontig$LabelChannel[y] == cntNick & 
                        (datcontig$LabelChannel[x] == cntNick 
                        & datcontig$contigPosition[x]>0)){
                        
                            g1<-grep(cntNick, as.character(datcontig$LabelChannel[1:z]))
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                                ji= length(g1)
                                while(ji > 0){
                                    z=g1[ji]
                                    #print (paste(z,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(contposn[z] > 0 
                                        & (datcontig$contigPosition[y]== 0 
                                        | is.nan(datcontig$contigPosition[y]))){
                                         datcontig$contigPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(datcontig$Position[x]) - as.numeric(posn[z])))* (datcontig$contigPosition[x] - contposn[z])) + contposn[z]
                                        ji= -1
                                    }
                                    else{
                                        ji=ji-1
                                        next
                                    }
                                }
                                if(datcontig$contigPosition[y]== 0 
                                | is.nan(datcontig$contigPosition[y])){
                                    if((orientMolecule == "+")){
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                           if(datcontig$contigPosition[y] > datcontig$contigPosition[x]
                                               | datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                      } else{
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                               if(datcontig$contigPosition[y] < 0){
                                                  datcontig$contigPosition[y] = 0
                                              }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    
                                       }
                                }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                            'if(orientMolecule == "+"){
                                #datcontig$contigRight[y]=as.numeric(datcontig$contigPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                
                                
                                
                                if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                    datcontig$contigPosition[y] <- datcontig$contigLeft[y]
                                }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigPosition[y] <- datcontig$contigRight[y]
                                }
                                else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                    datcontig$contigRight[y] = datcontig$contigRight[y]
                                }
                                 else{
                                    datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                                }
                                } else{
                                    datcontig$contigRight[y] = as.numeric(datcontig$contigPosition[x]) + (as.numeric(datcontig$Position[y])-as.numeric(datcontig$Position[x]))
                                    g1<-grep("1",as.character(datcontig$LabelChannel[1:z]))
                                    ji= length(g1)
                                    posn <- as.numeric(datcontig$Position[1:z])
                                    contposn <- as.numeric(datcontig$contigPosition[1:z])
                                    while(ji > 0){
                                        z=g1[ji]
                                        #print (paste(z,":",y,":",cmapID[ti],":",ji,sep=""))
                                        if(contposn[z]>0  & datcontig$contigPosition[y] == 0){
                                            if(orientMolecule == "+"){
                                                datcontig$contigLeft[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y])-as.numeric(posn[z]))
                                                datcontig$contigPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                             } else{
                                                 datcontig$contigLeft[y] = as.numeric(contposn[z]) - (as.numeric(posn[z])-as.numeric(datcontig$Position[y]))
                                                datcontig$contigPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                            }
                                            ji = -1
                                        }
                                        else{
                                            ji=ji - 1
                                            next
                                        }
                                    }
                                    if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                        datcontig$contigPosition[y] <- datcontig$contigLeft[y]
                                    }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                        datcontig$contigPosition[y] <- datcontig$contigRight[y]
                                    }
                                    else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                        datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                        datcontig$contigRight[y] = datcontig$contigRight[y]
                                    }
                                     else{
                                        datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                                    }
                                }       
                    '
                    }
                    else if(datcontig$LabelChannel[y] == cntNick & 
                        ((datcontig$LabelChannel[z] == cntMeth 
                        | datcontig$LabelChannel[x]== cntMeth) 
                        & (datcontig$contigPosition[z]==0
                        |datcontig$contigPosition[x]==0))){
                        g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                        g2 <- grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                        if((datcontig$contigPosition[y]== 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g1) >=1 & length(g2) >=1)){
                            ji=length(g1)
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                            while(ji > 0){
                                z=g1[ji]
                                #print(ji)
                                #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[z] > 0 
                                    & (datcontig$contigPosition[y]== 0 
                                    | is.nan(datcontig$contigPosition[y]))){
                                    
                                        #datcontig$contigLeft[y]=as.numeric(contposn[z]) + (as.numeric(posn[z]) - as.numeric(datcontig$Position[y]))
                                        ###Calculating x
                                        g11<-grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                                        posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                        contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                        ki=1
                                        while(ki <=length(g11)){
                                            x1=g11[ki]
                                            ##print (paste(x1,":",y,":",cmapID[ti],":",ki,sep=""))
                                            if(contposn1[x1] > 0 
                                                & (datcontig$contigPosition[y]== 0 
                                                | is.nan(datcontig$contigPosition[y]))){
                                                datcontig$contigPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(posn1[x1]) - as.numeric(posn[z])))* (contposn1[x1] - contposn[z]))+ contposn[z]
                                                ki=length(g11)+1
                                                ji =  -1
                                            }
                                            else{
                                                ki=ki+1
                                                next
                                            }
                                        }
                                    if(datcontig$contigPosition[y]== 0 
                                        | is.nan(datcontig$contigPosition[y])){
                                        if((orientMolecule == "+")){
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$contigPosition[y] < datcontig$contigPosition[z]| datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        } else{
                                           datcontig$contigPosition[y] = as.numeric(datcontig$contigPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                               if(datcontig$contigPosition[y] < 0){
                                                  datcontig$contigPosition[y] = 0
                                              }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        
                                        }
                                    }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}    
                                     
                                    
                                    ji =  ji - 1
                                    #print(paste("1:",ji))
                                }
                                else{
                                    datcontig$contigPosition[y]=0
                                    ji = ji-1
                                   # print(paste("2:",ji))
                                }
                            }
                            if(datcontig$contigPosition[y]== 0 
                                | is.nan(datcontig$contigPosition[y])){
                                ji=length(g1)
                                g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                                posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                                ji=1
                                while(ji <=length(g1)){
                                    x=g1[ji]
                                    #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(y < x & datcontig$contigPosition[x] > 0){
                                        if(orientMolecule == "+"){
                                            datcontig$contigPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                            if(datcontig$contigPosition[y] > contposn1[x] 
                                                | datcontig$contigPosition[y] < 0){
                                                datcontig$contigPosition[y] = 0
                                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    } else{
                                            datcontig$contigPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                            if(datcontig$contigPosition[y] < 0){
                                                datcontig$contigPosition[y] = 0
                                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                        }
                                        
                                        ji = ji+1
                                    }
                                    else{
                                        ji=ji+1
                                        next
                                    }
                                }
                            }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}    
                            
                        } else if((datcontig$contigPosition[y] == 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g2) >=1 & length(g1) == 0)){
                                ji=length(g1)
                                g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                                posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                contposn1 <- as.numeric(datcontig$contigPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y < x & datcontig$contigPosition[x] > 0){
                                    if(orientMolecule == "+"){
                                        datcontig$contigPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$contigPosition[y] > contposn1[x] 
                                            | datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                } else{
                                        datcontig$contigPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                    }
                                    
                                    ji = ji+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }            
            
                        } else if((datcontig$contigPosition[y]== 0 
                            | is.nan(datcontig$contigPosition[y]))
                            & (length(g2) == 0 & length(g1) >= 1)){
                            g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                            ji=length(g1)
                            posn <- as.numeric(datcontig$Position[1:z])
                            contposn <- as.numeric(datcontig$contigPosition[1:z])
                            while(ji > 0){
                                z=g1[ji]
                                #print(ji)
                                #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y > z & (contposn[z] > 0 & !is.na(contposn[z]))){
                                    if((orientMolecule == "+")){
                                        datcontig$contigPosition[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                        if(datcontig$contigPosition[y] < contposn[z]| datcontig$contigPosition[y] < 0){
                                            datcontig$contigPosition[y] = 0
                                        }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                   } else{
                                        datcontig$contigPosition[y] = as.numeric(contposn[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                            if(datcontig$contigPosition[y] < 0){
                                               datcontig$contigPosition[y] = 0
                                           }else{datcontig$contigPosition[y] = datcontig$contigPosition[y]}
                                
                                    }
                                       ji= ji - 1
                                } else{
                                    ji = ji - 1
                                    next
                                }               
                            }
                        }else{
                            #print("Locations do not match!!!!!")
                        }    
                    }else if (datcontig$LabelChannel[y] == cntNick & datcontig$contigPosition[y]>0){
                        datcontig$contigPosition[y]=datcontig$contigPosition[y]
                    }else{
                        datcontig$contigPosition[y]=0
                    
                    }
                    
                }
                if(datcontig$contigPosition[y] == 0 
                    & datcontig$LabelChannel[y] == 2){
                #print(paste("CMapID", cmapID[ti], "y:", y, sep = ""))
                'print(y)
                print(datcontig$contigPosition[y])'
                }else{
                   
                }
                 
            }
            
                
                dataContig<-rbind(dataContig,datcontig,row.names = NULL)
            
        }
	#stopImplicitCluster(cl)	
        'if(outputType == "text"){
            fname=paste("Contig_",contigID[ii],"MethPositions.txt",sep="")
            #rownames(datacontig)<-NULL
            print(file.path(hmdir,fname,fsep=""))
            write.table(dataContig,file.path(hmdir,fname,fsep="/"),row.names=FALSE)
        }
        else{return(dataContig)}'
		return(dataContig)
}
#' Liftover of molecule methylation and nick  contig position to reference coordinates
#'
#' @param molcmap character. molecule cmap data.
#' @param outputType character. Output method. dataframe or text.
#' @param cntNick integer. Nick label value.
#' @param cntMeth integer. Meth label value.
#' @return dataframe containing the methyl label positions mapped to the reference.
#' @examples
#' xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
#' cmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
#' modcmap <- readingCmap(cmap)
#' modxmap <- readingXmap(xmap)
#' modMolcmap <- modmolcmap(molcmap = modcmap,
#'   xmapdata = modxmap,
#'   cntNick = 1,
#'   cntMeth=2)
#' refcontigXmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/ContigRef.xmap", package="methometR")
#' refCmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/hg19ref_r.cmap", package="methometR")
#' Contigqcmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/SampContig_q.cmap", package="methometR")
#' refxmapdat <- readingXmap(refcontigXmap)
#' refcmapdat <- readingCmap(refCmap)
#' contigcmapdat <- readingCmap(Contigqcmap)
#' nickRef<-nickReference(refxmap = refxmapdat, refcmap = refcmapdat, 
#'     contigcmap = contigcmapdat, contigID = 6701, 
#'    returnMethod =c("dataFrame"),  chrom = 4)
#' xmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContigMolecule.xmap", package="methometR")
#' molcmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampMolecule_q.cmap", package="methometR")
#' contigcmap <- system.file("extdata",  "Cmapdir/output/contigs/exp_refineFinal1/merged_smaps/SampContig_r.cmap", package="methometR")
#' contigID = 6701; cntNick = 1; cntMeth=2
#' contigmodcmap <- readingXmap(contigcmap)
#' modmolcmapNick <-molNickLoc(nickrefLoc = nickRef,
#'    molcmap = modMolcmap,
#'    xmapdata = modxmap, 
#'    contigID = contigID, cntNick = cntNick, 
#'    cntMeth = cntMeth)
#' MethylRefLoc<-molMethylRefLoc(molcmap = modmolcmapNick,
#'    outputType = c("dataframe"),
#'    cntNick = cntNick, 
#'    cntMeth = cntMeth)
#' @import utils
#' @export
molMethylRefLoc <- function(molcmap, outputType = c("dataframe", "text"), cntMeth = 1, cntNick = 2){
    print("Searching for the methylation label locations in the reference")
    molcmap = molcmap
	#molcmap = molcmapFinal
    dataContig <- data.frame()
    datFinalmol = molcmap
    cmapID<-as.character(unique(datFinalmol$CMapId))
    #print(paste("cmapID::",cmapID))
    ###Adding the probable location of Methylation
    #datcontig<-data.frame()
	'cores=detectCores()
    cl <- makeCluster(5) #not to overload your computer
	print(cl)
    registerDoParallel(cl)'
	#dataContig <- foreach(ti=1:length(cmapID), .combine=rbind) %dopar% {
    for ( ti in 1:length(cmapID)){
        #print(paste("ti", ti))
        datcontig<-data.frame()
		#datcontig<-datFinalmol[which(as.numeric(datFinalmol$CMapId)== cmapID[ti]),]
        datcontig<-datFinalmol[which(as.numeric(datFinalmol$CMapId)== as.numeric(cmapID[ti])),]
        orientMolecule <- unique(as.character(datcontig$OrientationMolecule))
        #print(orientMolecule)
        datcontig[is.na(datcontig)]<-0
        
        for( y in 1:nrow(datcontig)){
            z=y-1
            x=y+1
            #print(paste("y", y))
            if(y == 1){
                if(datcontig$LabelChannel[y] == cntMeth 
                    & (datcontig$LabelChannel[x] == cntNick 
                    & datcontig$chromPosition[x]>0)){
                    if(orientMolecule == "+"){
                        datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                        if(datcontig$chromPosition[y] > datcontig$chromPosition[x]| datcontig$chromPosition[y] < 0){
                            datcontig$chromPosition[y] = 0
                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                    }else{
                        datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                        if(datcontig$chromPosition[y] < 0){
                            datcontig$chromPosition[y] = 0
                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                    
                    }
                }
                else if((datcontig$LabelChannel[y]== cntMeth 
                    & datcontig$chromPosition[y] == 0) 
                    & (datcontig$LabelChannel[x] == cntNick 
                    & datcontig$chromPosition[x] == 0)){
                    g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                    ji=1
                    while(ji <=length(g1)){
                        x=g1[ji]
                        #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                        if(y < x & datcontig$chromPosition[x] > 0){
                            if(orientMolecule == "+"){
                                datcontig$chromPosition[y]=as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                if(datcontig$chromPosition[y] > datcontig$chromPosition[x]| datcontig$chromPosition[y] < 0){
                                    datcontig$chromPosition[y] = 0
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                        } else{
                                datcontig$chromPosition[y]=as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                if(datcontig$chromPosition[y] < 0){
                                    datcontig$chromPosition[y] = 0
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            }
                            
                            ji=length(g1)+1
                        }
                        else{
                            ji=ji+1
                            next
                        }
                    }            
                
                }
                else if(datcontig$LabelChannel[y]== cntMeth 
                    & (datcontig$LabelChannel[x] == cntMeth 
                    & datcontig$chromPosition[x]==0)){
                    g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                    ji=1
                    while(ji <=length(g1)){
                        x=g1[ji]
                        #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                        if(y < x & datcontig$chromPosition[x] > 0){
                            if((orientMolecule == "+")){
                                datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                if(datcontig$chromPosition[y] > datcontig$chromPosition[x]| datcontig$chromPosition[y] < 0){
                                    datcontig$chromPosition[y] = 0
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            } else{
                                datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                if(datcontig$chromPosition[y] < 0){
                                    datcontig$chromPosition[y] = 0
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                            }
                        ji=length(g1)+1
                        }
                        else{
                            ji=ji+1
                            next
                        }
                    }
                }
                else if(datcontig$LabelChannel[y] == cntNick 
                    & (datcontig$LabelChannel[x] == cntNick 
                    & datcontig$chromPosition[x]>0)){
                    if(orientMolecule == "+"){
                        datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                        if(datcontig$chromPosition[y] > datcontig$chromPosition[x]| datcontig$chromPosition[y] < 0){
                            datcontig$chromPosition[y] = 0
                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                    }else{
                        datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                        if(datcontig$chromPosition[y] < 0){
                            datcontig$chromPosition[y] = 0
                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                    
                    }
                }
                else if((datcontig$LabelChannel[y]== cntNick 
                    & datcontig$chromPosition[y] == 0) 
                    & (datcontig$LabelChannel[x] == cntNick 
                    & datcontig$chromPosition[x] == 0)){
                    g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                    ji=1
                    while(ji <=length(g1)){
                        x=g1[ji]
                        #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                        if(y < x & datcontig$chromPosition[x] > 0){
                            if(orientMolecule == "+"){
                                datcontig$chromPosition[y]=as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                if(datcontig$chromPosition[y] > datcontig$chromPosition[x]| datcontig$chromPosition[y] < 0){
                                    datcontig$chromPosition[y] = 0
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                        } else{
                                datcontig$chromPosition[y]=as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                                if(datcontig$chromPosition[y] < 0){
                                    datcontig$chromPosition[y] = 0
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            }
                            
                            ji=length(g1)+1
                        }
                        else{
                            ji=ji+1
                            next
                        }
                    }            
                
                }
                else if(datcontig$LabelChannel[y] == cntNick 
                    & (datcontig$LabelChannel[x] == cntMeth 
                    & datcontig$chromPosition[x] == 0)){
                    g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                    ji=1
                    while(ji <=length(g1)){
                        x=g1[ji]
                        #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                        if(y < x & datcontig$chromPosition[x] > 0){
                            if((orientMolecule == "+")){
                                datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                if(datcontig$chromPosition[y] > datcontig$chromPosition[x]| datcontig$chromPosition[y] < 0){
                                    datcontig$chromPosition[y] = 0
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            } else{
                                datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                if(datcontig$chromPosition[y] < 0){
                                    datcontig$chromPosition[y] = 0
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                            }
                        ji=length(g1)+1
                        }
                        else{
                            ji=ji+1
                            next
                        }
                    }
                }
                else if (datcontig$LabelChannel[y] == cntNick & datcontig$chromPosition[y]>0){
                    datcontig$chromPosition[y] = datcontig$chromPosition[y]
                }
                else {
                
                    datcontig$chromPosition[y]=0
                }
            }
            else{
                if(datcontig$LabelChannel[y] == cntMeth 
                    & (datcontig$LabelChannel[z] == cntNick 
                    & datcontig$chromPosition[z]>0)){
                
                    'if(orientMolecule == "+"){'
                        'datcontig$contigLeft[y] <- as.numeric(datcontig$chromPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))'
                        g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:nrow(datcontig)]))
                        posn <- as.numeric(datcontig$Position[y+1:nrow(datcontig)])
                        contposn <- as.numeric(datcontig$chromPosition[y+1:nrow(datcontig)])
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(contposn[x] > 0  
                                & (datcontig$chromPosition[y]== 0 
                                | is.nan(datcontig$chromPosition[y]))){
                                datcontig$chromPosition[y] <- (
                    (
            (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))/(as.numeric(posn[x]) - as.numeric(datcontig$Position[z])))* (contposn[x] - datcontig$chromPosition[z]))+ datcontig$chromPosition[z]
                            ji=length(g1)+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }
                        if(datcontig$chromPosition[y]== 0 
                            | is.nan(datcontig$chromPosition[y])){
                            if((orientMolecule == "+")){
                                   datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                   if(datcontig$chromPosition[y] < datcontig$chromPosition[z]| datcontig$chromPosition[y] < 0){
                                       datcontig$chromPosition[y] = 0
                                   }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                              } else{
                                   datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                       if(datcontig$chromPosition[y] < 0){
                                          datcontig$chromPosition[y] = 0
                                      }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                               }
                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                        'if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                            datcontig$chromPosition[y] <- datcontig$contigLeft[y]
                        }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                            datcontig$chromPosition[y] <- datcontig$contigRight[y]
                        }
                        else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                            datcontig$contigLeft[y] = datcontig$contigLeft[y]
                            datcontig$contigRight[y] = datcontig$contigRight[y]
                        }
                        else{
                            datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                        }'
                    
                           '} else{
                            datcontig$contigLeft[y] = as.numeric(datcontig$chromPosition[z]) - (as.numeric(datcontig$Position[z])-as.numeric(datcontig$Position[y]))
                            g1<-grep(cntNick, as.character(datcontig$LabelChannel[(y+1):length(datcontig$LabelChannel)]))
                            posn <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[x]>0  & datcontig$chromPosition[y] == 0){
                                    if(orientMolecule == "+"){
                                        datcontig$contigRight[y] = as.numeric(contposn[x]) - (as.numeric(posn[x])-as.numeric(datcontig$Position[y]))
                                        datcontig$chromPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                     } else{
                                         datcontig$contigRight[y] = as.numeric(contposn[x]) + (as.numeric(datcontig$Position[y])-as.numeric(posn[x]))
                                        datcontig$chromPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                    }
                                ji=length(g1)+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                } 
                            }
                            
                            if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                datcontig$chromPosition[y] <- datcontig$contigLeft[y]
                            }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                datcontig$chromPosition[y] <- datcontig$contigRight[y]
                            }
                            else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                datcontig$contigRight[y] = datcontig$contigRight[y]
                            }
                             else{
                                datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                            }
                        }'
                }
                else if(datcontig$LabelChannel[y] == cntMeth & 
                    (datcontig$LabelChannel[x] == cntNick 
                    & datcontig$chromPosition[x]>0)){
                    
                        g1<-grep(cntNick, as.character(datcontig$LabelChannel[1:z]))
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                            ji= length(g1)
                            while(ji > 0){
                                z=g1[ji]
                                #print (paste(z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[z] > 0 
                                    & (datcontig$chromPosition[y]== 0 
                                    | is.nan(datcontig$chromPosition[y]))){
                                     datcontig$chromPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(datcontig$Position[x]) - as.numeric(posn[z])))* (datcontig$chromPosition[x] - contposn[z])) + contposn[z]
                                    ji= -1
                                }
                                else{
                                    ji=ji-1
                                    next
                                }
                            }
                            if(datcontig$chromPosition[y]== 0 
                            | is.nan(datcontig$chromPosition[y])){
                                if((orientMolecule == "+")){
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                       if(datcontig$chromPosition[y] > datcontig$chromPosition[x]
                                           | datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                  } else{
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                           if(datcontig$chromPosition[y] < 0){
                                              datcontig$chromPosition[y] = 0
                                          }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                
                                   }
                            }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                        'if(orientMolecule == "+"){
                            #datcontig$contigRight[y]=as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                            
                            
                            
                            if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                datcontig$chromPosition[y] <- datcontig$contigLeft[y]
                            }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                datcontig$chromPosition[y] <- datcontig$contigRight[y]
                            }
                            else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                datcontig$contigRight[y] = datcontig$contigRight[y]
                            }
                             else{
                                datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                            }
                            } else{
                                datcontig$contigRight[y] = as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[y])-as.numeric(datcontig$Position[x]))
                                g1<-grep("1",as.character(datcontig$LabelChannel[1:z]))
                                ji= length(g1)
                                posn <- as.numeric(datcontig$Position[1:z])
                                contposn <- as.numeric(datcontig$chromPosition[1:z])
                                while(ji > 0){
                                    z=g1[ji]
                                    #print (paste(z,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(contposn[z]>0  & datcontig$chromPosition[y] == 0){
                                        if(orientMolecule == "+"){
                                            datcontig$contigLeft[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y])-as.numeric(posn[z]))
                                            datcontig$chromPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                         } else{
                                             datcontig$contigLeft[y] = as.numeric(contposn[z]) - (as.numeric(posn[z])-as.numeric(datcontig$Position[y]))
                                            datcontig$chromPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                        }
                                        ji = -1
                                    }
                                    else{
                                        ji=ji - 1
                                        next
                                    }
                                }
                                if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                    datcontig$chromPosition[y] <- datcontig$contigLeft[y]
                                }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                    datcontig$chromPosition[y] <- datcontig$contigRight[y]
                                }
                                else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                    datcontig$contigRight[y] = datcontig$contigRight[y]
                                }
                                 else{
                                    datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                                }
                            }       
                '
                }
                else if(datcontig$LabelChannel[y] == cntMeth 
                    & ((datcontig$LabelChannel[z] == cntNick 
                    | datcontig$LabelChannel[x] == cntNick) 
                    & (datcontig$chromPosition[z] == 0
                    |datcontig$chromPosition[x]==0))){
                    g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                    g2 <- grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                    if((datcontig$chromPosition[y]== 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g1) >=1 & length(g2) >=1)){
                        ji=length(g1)
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                        while(ji > 0){
                            z=g1[ji]
                            'print(ji)
                            print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))'
                            if(contposn[z] > 0 
                                & (datcontig$chromPosition[y]== 0 
                                | is.nan(datcontig$chromPosition[y]))){
                                
                                    #datcontig$contigLeft[y]=as.numeric(contposn[z]) + (as.numeric(posn[z]) - as.numeric(datcontig$Position[y]))
                                    ###Calculating x
                                    g11<-grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                                    posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                    contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                                    ki=1
                                    while(ki <=length(g11)){
                                        x1=g11[ki]
                                        ##print (paste(x1,":",y,":",cmapID[ti],":",ki,sep=""))
                                        if(contposn1[x1] > 0 
                                            & (datcontig$chromPosition[y]== 0 
                                            | is.nan(datcontig$chromPosition[y]))){
                                            datcontig$chromPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(posn1[x1]) - as.numeric(posn[z])))* (contposn1[x1] - contposn[z]))+ contposn[z]
                                            ki=length(g11)+1
                                            ji =  -1
                                        }
                                        else{
                                            ki=ki+1
                                            next
                                        }
                                    }
                                if(datcontig$chromPosition[y]== 0 
                                    | is.nan(datcontig$chromPosition[y])){
                                    if((orientMolecule == "+")){
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                       if(datcontig$chromPosition[y] < datcontig$chromPosition[z]| datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    } else{
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$chromPosition[y] < 0){
                                              datcontig$chromPosition[y] = 0
                                          }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    
                                    }
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}    
                                 
                                
                                ji =  ji - 1
                                #print(paste("1:",ji))
                            }
                            else{
                                datcontig$chromPosition[y]=0
                                ji = ji-1
                                #print(paste("2:",ji))
                            }
                        }
                        if(datcontig$chromPosition[y]== 0 
                            | is.nan(datcontig$chromPosition[y])){
                            ji=length(g1)
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                            posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y < x & datcontig$chromPosition[x] > 0){
                                    if(orientMolecule == "+"){
                                        datcontig$chromPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$chromPosition[y] > contposn1[x] 
                                            | datcontig$chromPosition[y] < 0){
                                            datcontig$chromPosition[y] = 0
                                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                } else{
                                        datcontig$chromPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$chromPosition[y] < 0){
                                            datcontig$chromPosition[y] = 0
                                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    }
                                    
                                    ji = ji+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }
                        }else{
                            datcontig$chromPosition[y] = datcontig$chromPosition[y]
                        }    
                        
                    } else if((datcontig$chromPosition[y] == 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g2) >=1 & length(g1) == 0)){
                            ji=length(g1)
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                            posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(y < x & datcontig$chromPosition[x] > 0){
                                if(orientMolecule == "+"){
                                    datcontig$chromPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$chromPosition[y] > contposn1[x] 
                                        | datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            } else{
                                    datcontig$chromPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                }
                                
                                ji = ji+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }            
        
                    } else if((datcontig$chromPosition[y]== 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g2) == 0 & length(g1) >= 1)){
                        g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                        ji=length(g1)
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                        while(ji > 0){
                            z=g1[ji]
                            'print(ji)
                            print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))'
                            if(y > z & (contposn[z] > 0 & !is.na(contposn[z]))){
                                if((orientMolecule == "+")){
                                    datcontig$chromPosition[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                    if(datcontig$chromPosition[y] < contposn[z]| datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                               } else{
                                    datcontig$chromPosition[y] = as.numeric(contposn[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                        if(datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                                }
                                   ji= ji - 1
                            } else{
                                ji = ji - 1
                                next
                            }               
                        }
                    }else{
                        #print("Locations do not match!!!!!")
                    }    
                }
                else if(datcontig$LabelChannel[y] == cntMeth & 
                    ((datcontig$LabelChannel[z] == cntMeth 
                    | datcontig$LabelChannel[x]== cntMeth) 
                    & (datcontig$chromPosition[z]==0
                    |datcontig$chromPosition[x]==0))){g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                    g2 <- grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                    if((datcontig$chromPosition[y]== 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g1) >=1 & length(g2) >=1)){
                        ji=length(g1)
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                        while(ji > 0){
                            z=g1[ji]
                            'print(ji)
                            print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))'
                            if(contposn[z] > 0 
                                & (datcontig$chromPosition[y]== 0 
                                | is.nan(datcontig$chromPosition[y]))){
                                
                                    #datcontig$contigLeft[y]=as.numeric(contposn[z]) + (as.numeric(posn[z]) - as.numeric(datcontig$Position[y]))
                                    ###Calculating x
                                    g11<-grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                                    posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                    contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                                    ki=1
                                    while(ki <=length(g11)){
                                        x1=g11[ki]
                                        ##print (paste(x1,":",y,":",cmapID[ti],":",ki,sep=""))
                                        if(contposn1[x1] > 0 
                                            & (datcontig$chromPosition[y]== 0 
                                            | is.nan(datcontig$chromPosition[y]))){
                                            datcontig$chromPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(posn1[x1]) - as.numeric(posn[z])))* (contposn1[x1] - contposn[z]))+ contposn[z]
                                            ki=length(g11)+1
                                            ji =  -1
                                        }
                                        else{
                                            ki=ki+1
                                            next
                                        }
                                    }
                                if(datcontig$chromPosition[y]== 0 
                                    | is.nan(datcontig$chromPosition[y])){
                                    if((orientMolecule == "+")){
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                       if(datcontig$chromPosition[y] < datcontig$chromPosition[z]| datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    } else{
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$chromPosition[y] < 0){
                                              datcontig$chromPosition[y] = 0
                                          }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    
                                    }
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}    
                                 
                                
                                ji =  ji - 1
                                #print(paste("1:",ji))
                            }
                            else{
                                datcontig$chromPosition[y]=0
                                ji = ji-1
                                #print(paste("2:",ji))
                            }
                        }
                        if(datcontig$chromPosition[y]== 0 
                            | is.nan(datcontig$chromPosition[y])){
                            ji=length(g1)
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                            posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y < x & datcontig$chromPosition[x] > 0){
                                    if(orientMolecule == "+"){
                                        datcontig$chromPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$chromPosition[y] > contposn1[x] 
                                            | datcontig$chromPosition[y] < 0){
                                            datcontig$chromPosition[y] = 0
                                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                } else{
                                        datcontig$chromPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$chromPosition[y] < 0){
                                            datcontig$chromPosition[y] = 0
                                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    }
                                    
                                    ji = ji+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }
                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}    
                        
                    } else if((datcontig$chromPosition[y] == 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g2) >=1 & length(g1) == 0)){
                            ji=length(g1)
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                            posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(y < x & datcontig$chromPosition[x] > 0){
                                if(orientMolecule == "+"){
                                    datcontig$chromPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$chromPosition[y] > contposn1[x] 
                                        | datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            } else{
                                    datcontig$chromPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                }
                                
                                ji = ji+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }            
        
                    } else if((datcontig$chromPosition[y]== 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g2) == 0 & length(g1) >= 1)){
                        g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                        ji=length(g1)
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                        while(ji > 0){
                            z=g1[ji]
                            'print(ji)
                            print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))'
                            if(y > z & (contposn[z] > 0 & !is.na(contposn[z]))){
                                if((orientMolecule == "+")){
                                    datcontig$chromPosition[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                    if(datcontig$chromPosition[y] < contposn[z]| datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                               } else{
                                    datcontig$chromPosition[y] = as.numeric(contposn[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                        if(datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                                }
                                   ji= ji - 1
                            } else{
                                ji = ji - 1
                                next
                            }               
                        }
                    }else{
                        #print("Locations do not match!!!!!")
                    }    
                }
                else if(datcontig$LabelChannel[y] == cntNick 
                    & (datcontig$LabelChannel[z] == cntNick 
                    & datcontig$chromPosition[z]>0)){
                
                    'if(orientMolecule == "+"){'
                        'datcontig$contigLeft[y] <- as.numeric(datcontig$chromPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))'
                        g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:nrow(datcontig)]))
                        posn <- as.numeric(datcontig$Position[y+1:nrow(datcontig)])
                        contposn <- as.numeric(datcontig$chromPosition[y+1:nrow(datcontig)])
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(contposn[x] > 0  
                                & (datcontig$chromPosition[y]== 0 
                                | is.nan(datcontig$chromPosition[y]))){
                                datcontig$chromPosition[y] <- (
                    (
            (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))/(as.numeric(posn[x]) - as.numeric(datcontig$Position[z])))* (contposn[x] - datcontig$chromPosition[z]))+ datcontig$chromPosition[z]
                            ji=length(g1)+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }
                        if(datcontig$chromPosition[y]== 0 
                            | is.nan(datcontig$chromPosition[y])){
                            if((orientMolecule == "+")){
                                   datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                   if(datcontig$chromPosition[y] < datcontig$chromPosition[z]| datcontig$chromPosition[y] < 0){
                                       datcontig$chromPosition[y] = 0
                                   }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                              } else{
                                   datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                       if(datcontig$chromPosition[y] < 0){
                                          datcontig$chromPosition[y] = 0
                                      }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                               }
                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                        'if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                            datcontig$chromPosition[y] <- datcontig$contigLeft[y]
                        }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                            datcontig$chromPosition[y] <- datcontig$contigRight[y]
                        }
                        else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                            datcontig$contigLeft[y] = datcontig$contigLeft[y]
                            datcontig$contigRight[y] = datcontig$contigRight[y]
                        }
                        else{
                            datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                        }'
                    
                           '} else{
                            datcontig$contigLeft[y] = as.numeric(datcontig$chromPosition[z]) - (as.numeric(datcontig$Position[z])-as.numeric(datcontig$Position[y]))
                            g1<-grep(cntNick, as.character(datcontig$LabelChannel[(y+1):length(datcontig$LabelChannel)]))
                            posn <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[x]>0  & datcontig$chromPosition[y] == 0){
                                    if(orientMolecule == "+"){
                                        datcontig$contigRight[y] = as.numeric(contposn[x]) - (as.numeric(posn[x])-as.numeric(datcontig$Position[y]))
                                        datcontig$chromPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                     } else{
                                         datcontig$contigRight[y] = as.numeric(contposn[x]) + (as.numeric(datcontig$Position[y])-as.numeric(posn[x]))
                                        datcontig$chromPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                    }
                                ji=length(g1)+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                } 
                            }
                            
                            if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                datcontig$chromPosition[y] <- datcontig$contigLeft[y]
                            }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                datcontig$chromPosition[y] <- datcontig$contigRight[y]
                            }
                            else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                datcontig$contigRight[y] = datcontig$contigRight[y]
                            }
                             else{
                                datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                            }
                        }'
                }
                else if(datcontig$LabelChannel[y] == cntNick & 
                    (datcontig$LabelChannel[x] == cntNick 
                    & datcontig$chromPosition[x]>0)){
                    
                        g1<-grep(cntNick, as.character(datcontig$LabelChannel[1:z]))
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                            ji= length(g1)
                            while(ji > 0){
                                z=g1[ji]
                                #print (paste(z,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(contposn[z] > 0 
                                    & (datcontig$chromPosition[y]== 0 
                                    | is.nan(datcontig$chromPosition[y]))){
                                     datcontig$chromPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(datcontig$Position[x]) - as.numeric(posn[z])))* (datcontig$chromPosition[x] - contposn[z])) + contposn[z]
                                    ji= -1
                                }
                                else{
                                    ji=ji-1
                                    next
                                }
                            }
                            if(datcontig$chromPosition[y]== 0 
                            | is.nan(datcontig$chromPosition[y])){
                                if((orientMolecule == "+")){
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                       if(datcontig$chromPosition[y] > datcontig$chromPosition[x]
                                           | datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                  } else{
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[x]) - as.numeric(datcontig$Position[y]))
                                           if(datcontig$chromPosition[y] < 0){
                                              datcontig$chromPosition[y] = 0
                                          }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                
                                   }
                            }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                        'if(orientMolecule == "+"){
                            #datcontig$contigRight[y]=as.numeric(datcontig$chromPosition[x]) - (as.numeric(datcontig$Position[x])-as.numeric(datcontig$Position[y]))
                            
                            
                            
                            if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                datcontig$chromPosition[y] <- datcontig$contigLeft[y]
                            }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                datcontig$chromPosition[y] <- datcontig$contigRight[y]
                            }
                            else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                datcontig$contigRight[y] = datcontig$contigRight[y]
                            }
                             else{
                                datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                            }
                            } else{
                                datcontig$contigRight[y] = as.numeric(datcontig$chromPosition[x]) + (as.numeric(datcontig$Position[y])-as.numeric(datcontig$Position[x]))
                                g1<-grep("1",as.character(datcontig$LabelChannel[1:z]))
                                ji= length(g1)
                                posn <- as.numeric(datcontig$Position[1:z])
                                contposn <- as.numeric(datcontig$chromPosition[1:z])
                                while(ji > 0){
                                    z=g1[ji]
                                    #print (paste(z,":",y,":",cmapID[ti],":",ji,sep=""))
                                    if(contposn[z]>0  & datcontig$chromPosition[y] == 0){
                                        if(orientMolecule == "+"){
                                            datcontig$contigLeft[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y])-as.numeric(posn[z]))
                                            datcontig$chromPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                         } else{
                                             datcontig$contigLeft[y] = as.numeric(contposn[z]) - (as.numeric(posn[z])-as.numeric(datcontig$Position[y]))
                                            datcontig$chromPosition[y] <- (datcontig$contigRight[y] + datcontig$contigLeft[y])/2
                                        }
                                        ji = -1
                                    }
                                    else{
                                        ji=ji - 1
                                        next
                                    }
                                }
                                if(datcontig$contigRight[y] == 0 & datcontig$contigLeft[y] > 0){
                                    datcontig$chromPosition[y] <- datcontig$contigLeft[y]
                                }else if (datcontig$contigLeft[y] == 0 & datcontig$contigRight[y] > 0){
                                    datcontig$chromPosition[y] <- datcontig$contigRight[y]
                                }
                                else if (datcontig$contigLeft[y] > 0 & datcontig$contigRight[y] > 0){
                                    datcontig$contigLeft[y] = datcontig$contigLeft[y]
                                    datcontig$contigRight[y] = datcontig$contigRight[y]
                                }
                                 else{
                                    datcontig$contigLeft[y] = 0; datcontig$contigRight[y] = 0;
                                }
                            }       
                '
                }                
                else if(datcontig$LabelChannel[y] == cntNick 
                    & ((datcontig$LabelChannel[z] == cntNick 
                    | datcontig$LabelChannel[x] == cntNick) 
                    & (datcontig$chromPosition[z] == 0
                    |datcontig$chromPosition[x]==0))){
                    g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                    g2 <- grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                    if((datcontig$chromPosition[y]== 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g1) >=1 & length(g2) >=1)){
                        ji=length(g1)
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                        while(ji > 0){
                            z=g1[ji]
                            'print(ji)
                            print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))'
                            if(contposn[z] > 0 
                                & (datcontig$chromPosition[y]== 0 
                                | is.nan(datcontig$chromPosition[y]))){
                                
                                    #datcontig$contigLeft[y]=as.numeric(contposn[z]) + (as.numeric(posn[z]) - as.numeric(datcontig$Position[y]))
                                    ###Calculating x
                                    g11<-grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                                    posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                    contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                                    ki=1
                                    while(ki <=length(g11)){
                                        x1=g11[ki]
                                        ##print (paste(x1,":",y,":",cmapID[ti],":",ki,sep=""))
                                        if(contposn1[x1] > 0 
                                            & (datcontig$chromPosition[y]== 0 
                                            | is.nan(datcontig$chromPosition[y]))){
                                            datcontig$chromPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(posn1[x1]) - as.numeric(posn[z])))* (contposn1[x1] - contposn[z]))+ contposn[z]
                                            ki=length(g11)+1
                                            ji =  -1
                                        }
                                        else{
                                            ki=ki+1
                                            next
                                        }
                                    }
                                if(datcontig$chromPosition[y]== 0 
                                    | is.nan(datcontig$chromPosition[y])){
                                    if((orientMolecule == "+")){
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                       if(datcontig$chromPosition[y] < datcontig$chromPosition[z]| datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    } else{
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$chromPosition[y] < 0){
                                              datcontig$chromPosition[y] = 0
                                          }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    
                                    }
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}    
                                 
                                
                                ji =  ji - 1
                                #print(paste("1:",ji))
                            }
                            else{
                                datcontig$chromPosition[y]=0
                                ji = ji-1
                                #print(paste("2:",ji))
                            }
                        }
                        if(datcontig$chromPosition[y]== 0 
                            | is.nan(datcontig$chromPosition[y])){
                            ji=length(g1)
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                            posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y < x & datcontig$chromPosition[x] > 0){
                                    if(orientMolecule == "+"){
                                        datcontig$chromPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$chromPosition[y] > contposn1[x] 
                                            | datcontig$chromPosition[y] < 0){
                                            datcontig$chromPosition[y] = 0
                                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                } else{
                                        datcontig$chromPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$chromPosition[y] < 0){
                                            datcontig$chromPosition[y] = 0
                                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    }
                                    
                                    ji = ji+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }
                        }else{
                            datcontig$chromPosition[y] = datcontig$chromPosition[y]
                        }    
                        
                    } else if((datcontig$chromPosition[y] == 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g2) >=1 & length(g1) == 0)){
                            ji=length(g1)
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                            posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(y < x & datcontig$chromPosition[x] > 0){
                                if(orientMolecule == "+"){
                                    datcontig$chromPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$chromPosition[y] > contposn1[x] 
                                        | datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            } else{
                                    datcontig$chromPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                }
                                
                                ji = ji+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }            
        
                    } else if((datcontig$chromPosition[y]== 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g2) == 0 & length(g1) >= 1)){
                        g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                        ji=length(g1)
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                        while(ji > 0){
                            z=g1[ji]
                           ' print(ji)
                            print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))'
                            if(y > z & (contposn[z] > 0 & !is.na(contposn[z]))){
                                if((orientMolecule == "+")){
                                    datcontig$chromPosition[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                    if(datcontig$chromPosition[y] < contposn[z]| datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                               } else{
                                    datcontig$chromPosition[y] = as.numeric(contposn[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                        if(datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                                }
                                   ji= ji - 1
                            } else{
                                ji = ji - 1
                                next
                            }               
                        }
                    }else{
                        #print("Locations do not match!!!!!")
                    }    
                }
                else if(datcontig$LabelChannel[y] == cntNick & 
                    ((datcontig$LabelChannel[z] == cntMeth 
                    | datcontig$LabelChannel[x]== cntMeth) 
                    & (datcontig$chromPosition[z]==0
                    |datcontig$chromPosition[x]==0))){
                    g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                    g2 <- grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                    if((datcontig$chromPosition[y]== 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g1) >=1 & length(g2) >=1)){
                        ji=length(g1)
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                        while(ji > 0){
                            z=g1[ji]
                            #print(ji)
                            #print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(contposn[z] > 0 
                                & (datcontig$chromPosition[y]== 0 
                                | is.nan(datcontig$chromPosition[y]))){
                                
                                    #datcontig$contigLeft[y]=as.numeric(contposn[z]) + (as.numeric(posn[z]) - as.numeric(datcontig$Position[y]))
                                    ###Calculating x
                                    g11<-grep(cntNick,as.character(datcontig$LabelChannel[(y+1):nrow(datcontig)]))
                                    posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                                    contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                                    ki=1
                                    while(ki <=length(g11)){
                                        x1=g11[ki]
                                        ##print (paste(x1,":",y,":",cmapID[ti],":",ki,sep=""))
                                        if(contposn1[x1] > 0 
                                            & (datcontig$chromPosition[y]== 0 
                                            | is.nan(datcontig$chromPosition[y]))){
                                            datcontig$chromPosition[y] <- (((as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))/(as.numeric(posn1[x1]) - as.numeric(posn[z])))* (contposn1[x1] - contposn[z]))+ contposn[z]
                                            ki=length(g11)+1
                                            ji =  -1
                                        }
                                        else{
                                            ki=ki+1
                                            next
                                        }
                                    }
                                if(datcontig$chromPosition[y]== 0 
                                    | is.nan(datcontig$chromPosition[y])){
                                    if((orientMolecule == "+")){
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                       if(datcontig$chromPosition[y] < datcontig$chromPosition[z]| datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    } else{
                                       datcontig$chromPosition[y] = as.numeric(datcontig$chromPosition[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(datcontig$Position[z]))
                                           if(datcontig$chromPosition[y] < 0){
                                              datcontig$chromPosition[y] = 0
                                          }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    
                                    }
                                }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}    
                                 
                                
                                ji =  ji - 1
                                #print(paste("1:",ji))
                            }
                            else{
                                datcontig$chromPosition[y]=0
                                ji = ji-1
                                #print(paste("2:",ji))
                            }
                        }
                        if(datcontig$chromPosition[y]== 0 
                            | is.nan(datcontig$chromPosition[y])){
                            ji=length(g1)
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                            posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                            ji=1
                            while(ji <=length(g1)){
                                x=g1[ji]
                                #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                                if(y < x & datcontig$chromPosition[x] > 0){
                                    if(orientMolecule == "+"){
                                        datcontig$chromPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$chromPosition[y] > contposn1[x] 
                                            | datcontig$chromPosition[y] < 0){
                                            datcontig$chromPosition[y] = 0
                                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                } else{
                                        datcontig$chromPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                        if(datcontig$chromPosition[y] < 0){
                                            datcontig$chromPosition[y] = 0
                                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                    }
                                    
                                    ji = ji+1
                                }
                                else{
                                    ji=ji+1
                                    next
                                }
                            }
                        }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}    
                        
                    } else if((datcontig$chromPosition[y] == 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g2) >=1 & length(g1) == 0)){
                            ji=length(g1)
                            g1<-grep(cntNick,as.character(datcontig$LabelChannel[y+1:length(datcontig$LabelChannel)]))
                            posn1 <- as.numeric(datcontig$Position[(y+1):nrow(datcontig)])
                            contposn1 <- as.numeric(datcontig$chromPosition[(y+1):nrow(datcontig)])
                        ji=1
                        while(ji <=length(g1)){
                            x=g1[ji]
                            #print (paste(x,":",y,":",cmapID[ti],":",ji,sep=""))
                            if(y < x & datcontig$chromPosition[x] > 0){
                                if(orientMolecule == "+"){
                                    datcontig$chromPosition[y]=as.numeric(contposn1[x]) - (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$chromPosition[y] > contposn1[x] 
                                        | datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            } else{
                                    datcontig$chromPosition[y]=as.numeric(contposn1[x]) + (as.numeric(posn1[x])-as.numeric(datcontig$Position[y]))
                                    if(datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                                }
                                
                                ji = ji+1
                            }
                            else{
                                ji=ji+1
                                next
                            }
                        }            
        
                    } else if((datcontig$chromPosition[y]== 0 
                        | is.nan(datcontig$chromPosition[y]))
                        & (length(g2) == 0 & length(g1) >= 1)){
                        g1 <- grep(cntNick,as.character(datcontig$LabelChannel[1:z]))
                        ji=length(g1)
                        posn <- as.numeric(datcontig$Position[1:z])
                        contposn <- as.numeric(datcontig$chromPosition[1:z])
                        while(ji > 0){
                            z=g1[ji]
                            'print(ji)
                            print (paste("V:",z,":",y,":",cmapID[ti],":",ji,sep=""))'
                            if(y > z & (contposn[z] > 0 & !is.na(contposn[z]))){
                                if((orientMolecule == "+")){
                                    datcontig$chromPosition[y] = as.numeric(contposn[z]) + (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                    if(datcontig$chromPosition[y] < contposn[z]| datcontig$chromPosition[y] < 0){
                                        datcontig$chromPosition[y] = 0
                                    }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                               } else{
                                    datcontig$chromPosition[y] = as.numeric(contposn[z]) - (as.numeric(datcontig$Position[y]) - as.numeric(posn[z]))
                                        if(datcontig$chromPosition[y] < 0){
                                           datcontig$chromPosition[y] = 0
                                       }else{datcontig$chromPosition[y] = datcontig$chromPosition[y]}
                            
                                }
                                   ji= ji - 1
                            } else{
                                ji = ji - 1
                                next
                            }               
                        }
                    }else{
                        #print("Locations do not match!!!!!")
                    }    
                }
                else if (datcontig$LabelChannel[y] == cntNick & datcontig$chromPosition[y]>0){
                    datcontig$chromPosition[y]=datcontig$chromPosition[y]
                }
                else{
                    datcontig$chromPosition[y]=0
                
                }
                
            }
             
        }
            dataContig<-rbind(dataContig,datcontig,row.names = NULL)
            
        }
	#stopImplicitCluster(cl)	
        'if(outputType == "text"){
            fname=paste("Contig_",contigID[ii],"MethPositions.txt",sep="")
            #rownames(datacontig)<-NULL
            print(file.path(hmdir,fname,fsep=""))
            write.table(dataContig,file.path(hmdir,fname,fsep="/"),row.names=FALSE)
        }
        else{return(dataContig)}'
		return(dataContig)
}
#' Extracting contig information
#'
#' @param contigRefxmap character. Contig xmap data.
#' @param chrom character. Chromosome of interest.
#' @param startPos integer.Start position.
#' @param endPos integer. End position.
#' @return dataframe with data extracted from xmap.
#' @examples
#' refcontigXmap <- system.file("extdata", "Cmapdir/output/contigs/exp_refineFinal1_sv/merged_smaps/ContigRef.xmap", package="methometR")
#' modxmap <- readingXmap(refcontigXmap)
#' extractedInfo<-contigextract(contigRefxmap = modxmap, chrom = 4, startPos = 100000, endPos = 200000)
#' @import utils
#' @export
contigextract <- function(contigRefxmap, chrom, startPos, endPos){
    print("Extracting Contig information")
    if(startPos > 0 & endPos > 0 & chrom >0){
            dat1 <- contigRefxmap[which(contigRefxmap$RefStartPos <= startPos 
                & contigRefxmap$RefEndPos >= endPos 
                & contigRefxmap$RefContigID == chrom ),]
            
        }else if(startPos > 0 & endPos == 0 & chrom >0){
            dat1 <- contigRefxmap[which(contigRefxmap$RefStartPos <= startPos 
                & contigRefxmap$RefEndPos >= startPos 
                & contigRefxmap$RefContigID == chrom ),]
        }else if(startPos == 0 & endPos == 0 & chrom >0){
            dat1 <- contigRefxmap[which(contigRefxmap$RefContigID == chrom ),]
        }else{
            print("chromosome or start or end information required!!")
        }
    contigID <- dat1$QryContigID 
    chrom <- unique(dat1$RefContigID)
    dat <- data.frame(contigID = contigID, chrom = rep(chrom, length(contigID)))
    return(dat)
}
#' Main function to call the extraction, modification and quantification process.
#'
#' @param cmaploc_mol character. path to the cmap folder.
#' @param untar boolean. untar option. True or False/
#' @param cntNick integer. Nick label value.
#' @param cntMeth integer. Meth label value.
#' @param selectedRegion boolean. if TRUE, selected region is analyzed and 
#'     the user provides chromosome and start and end point. If FALSE, whole
#'     genome analysed.
#' @param chrom character. Chromosome number if selectedRegion is TRUE
#' @param startPos integer. Start of the query if selectedRegion is TRUE
#' @param endPos integer. End of the query if selectedRegion is TRUE
#' @param outputFormat Character. Choices are individualText, multipleText and
#'    dataframe.
#' @param path Character. Path of the folder where log files needs to be written 
#'     to.
#' @param logFilename character. Log file name.
#' @return dataframe or text files containing the modified files.
#' @examples
#'\dontrun{
#' molcmap <- system.file("extdata", "SampMolecule.cmap", package="methometR")
#' cmap <- system.file("extdata", "SampContig.cmap", package="methometR")
#' bed<-molMethylRefLoc(nickrefLoc,molcmap,xmapdata, contigID, cntNick , 
#'    cntMeth)
#'}
#' @import utils
#' @import logr
#' @export
extracting_modifying_methLabels <- function(cmaploc_mol, 
    untar = FALSE, cntNick, cntMeth, selectedRegion = TRUE, chrom, startPos, endPos,
    outputFormat = c("individualText", "multipleText", 
        "dataFrame"), path, logFilename){
	#logFilename = ""
	logFilename = logFilename
	ifelse(logFilename == "", 
	    log_open(file.path(path, "logmethometR.txt")),
		log_open(file.path(path, logFilename)))
	log_print("----methometr Package----")
	print("----Function to map nick and methylation label to reference----")
    log_print("----Function to map nick and methylation label to reference----")
	cntNick =cntNick; cntMeth = cntMeth;
    untar = untar; cmaploc_mol = cmaploc_mol;
    selectedRegion = selectedRegion;
    outputFormat = outputFormat; 
	path = path
	refcmap <- data.frame();
	querycmap <- data.frame();contigRefxmap <- data.frame()
	
    if(untar == TRUE){
        print("---Untaring and  reading molecule and contig info---")
		log_print("---Untaring and  reading molecule and contig info---")
        untaringdirectoy(cmaploc_mol)
        st <-  strsplit(cmaploc_mol, split = ".zip")
        cmapPath <- st[[1]][1]
        ##contig and reference maps
		print("Extracting file paths for required contig and reference files")
		contigRefpath <- file.path(cmapPath, "output/contigs/exp_refineFinal1_sv/merged_smaps/") 
        contRefXmappath <- list.files(contigRefpath, pattern = "*.xmap", full.names = TRUE)
        Refcmappath <- list.files(contigRefpath, pattern = "*_r.cmap", full.names = TRUE)
		Querycmappath <- list.files(contigRefpath, pattern = "*_q.cmap", full.names = TRUE)
		print("Reading contig, molecule and reference xmap")
        refcmap <- readingCmap(Refcmappath)
        querycmap <- readingCmap(Querycmappath)
        contigRefxmap <- readingXmap(xmap = contRefXmappath)
        chrom1 <-  unique(contigRefxmap$RefContigID)
		##contig and molecule maps
        contigMolPath <- file.path(cmapPath, "output/contigs/exp_refineFinal1/alignmol/merge/") 
        lmolecules <- list.files(cmaploc_mol, pattern = "*.xmap", full.names = FALSE)        
        st1 <- strsplit(lmolecules, split = "contig")
        st2 <- sapply(st1, "[[", 2)
        st3 <- strsplit(st2, split = "\\.")
        contigID <- sapply(st3, "[[", 1)
        lcontig <- list.files(cmaploc_mol, pattern = "*_r.cmap", full.names = TRUE)
        lxmap <- list.files(cmaploc_mol, pattern = "*.", full.names = TRUE)
    }
    else {
	    print("---Reading molecule and contig info---")
		log_print("---Reading molecule and contig info---")
	    ##contig and reference maps
		print("Extracting file paths for required contig and reference files")
        contigRefpath <- file.path(cmaploc_mol, "output/contigs/exp_refineFinal1_sv/merged_smaps/") 
        contRefXmappath <- list.files(contigRefpath, pattern = "*.xmap", full.names = TRUE)
        Refcmappath <- list.files(contigRefpath, pattern = "*_r.cmap", full.names = TRUE)
		Querycmappath <- list.files(contigRefpath, pattern = "*_q.cmap", full.names = TRUE)
        print("Reading contig, molecule and reference xmap")
		refcmap <- readingCmap(Refcmappath)
        querycmap <- readingCmap(Querycmappath)
        contigRefxmap <- readingXmap(xmap = contRefXmappath)
		chrom1 <-  unique(contigRefxmap$RefContigID)
        ##contig and molecule maps
        contigMolPath <- file.path(cmaploc_mol, "output/contigs/exp_refineFinal1/alignmol/merge") 
        lmolecules <- list.files(contigMolPath, pattern = "*.xmap", full.names = FALSE)
                
        st1 <- strsplit(lmolecules, split = "contig")
        st2 <- sapply(st1, "[[", 2)
        st3 <- strsplit(st2, split = "\\.")
        contigID <- sapply(st3, "[[", 1)
    }
    datmolInfo <- data.frame()
    if(selectedRegion == FALSE){
        modmol <- data.frame()
		print("Selected region is FALSE")
		'pp <- paste("Total number of contigs to be processed:", length(contigID))
		log_print(pp)
		print(paste("Total number of contigs to be processed:", length(contigID)))'
		chrom <- sort(chrom1)
		for(uu in seq_along(chrom)){
		    datcont <- contigRefxmap[which(contigRefxmap$RefContigID == chrom[uu]),]
			contigID <- datcont$QryContigID
			pa <- paste("---Processing Chromosome", 
			    paste("chr_", chrom[uu], sep = ""), 
				"--", sep = "")
			print(pa)
			log_print(pa)
			pa1 <- paste("---Number of contigs for chromosme ", 
			    paste("chr_", chrom[uu], sep = ""), ":", 
				length(contigID), "---", sep = "")
			
			#pp1 <- paste("---Number of contigs for chromosme1:", length(contigID))
			print(pa1)
			log_print(pa1)
			if(dir.exists(file.path(path, paste("chr_", chrom[uu], sep = "")))
			    & length(list.files(path = file.path(path, paste("chr_", chrom[uu], sep = "")))) > 0){
			    pa7 <- paste("chr_", chrom[uu],": already processed", sep = "")
				log_print(pa7)
				next;
			}else if(!dir.exists(file.path(path, paste("chr_", chrom[uu], sep = "")))){
			    pp2 <- paste("---Creating directory for Chromosome", chrom[uu], "--", sep = "")
			    print(paste("---Creating directory for Chromosome", chrom[uu], "--", sep = ""))
			    log_print(pp2)
			    dir.create(file.path(path, paste("chr_", chrom[uu], sep = "")))
				outpath <- file.path(path, paste("chr_", chrom[uu], sep = ""))
			}else{
			    outpath <- file.path(path, paste("chr_", chrom[uu], sep = ""))
			}
            
            if(length(contigID) > 1){
                
                for (ii in seq_along(contigID)){
				    contigRefxmap = contigRefxmap;refcmap = refcmap
		    	    #for (ii in 803:length(contigID)){
		    	    print(paste("Contig being analyzed:", contigID[ii], sep = ""))
					pp5 <- paste("Contig being analyzed:", contigID[ii], sep = "")
					log_print(pp5)
                    print("Reading contig, molecule and reference xmap")
					contigmap <- list.files (contigMolPath, pattern = paste("exp_refineFinal1_contig",contigID[ii],"_r.cmap",sep = ""), full.names = TRUE)
                    molmap <- list.files (contigMolPath, pattern = paste("exp_refineFinal1_contig",contigID[ii],"_q.cmap",sep = ""), full.names = TRUE)
                    xmap <- list.files (contigMolPath, pattern = paste("exp_refineFinal1_contig",contigID[ii],".xmap",sep = ""), full.names = TRUE)
                    #contigdat <- readingCmap(contigmap)
                    if(molmap == ""){
                        print(paste("exp_refineFinal1_contig",contigID[ii],"_q.cmap", "is discarded from analysis!!!"))
                        pp4 <- paste("exp_refineFinal1_contig",contigID[ii],"_q.cmap", "is discarded from analysis!!!")
						log_print(pp4)
						next;
                    }else{
                        moldat <- readingCmap(molmap)
                    }
                    dat1_t <- contigRefxmap[which(contigRefxmap$QryContigID == contigID[ii]),]
					if((nrow(dat1_t) > 1) 
					    & (length(unique(dat1_t$RefContigID)) > 1)){
						message("Work in progress translocation function")
						next;
					}else{
					    dat1 <- contigRefxmap[which(contigRefxmap$QryContigID == contigID[ii]
						    & contigRefxmap$RefContigID == chrom[uu]),]
					}
                    xmapdat <- readingXmap(xmap)
                    if(all(is.na((moldat[1:3,]))) == TRUE 
		    		    | all(is.na((xmapdat[1:3,]))) == TRUE  
		    			| all(is.na((querycmap[1:3,]))) == TRUE){
						log_print("Either xmap,molecule cmap or reference cmap file empty !!!!")
		    			print("Either xmap,molecule cmap or reference cmap file empty !!!!")
		    			next;
		    		}else{moldat = moldat; xmapdat = xmapdat}
                    if(nrow(dat1) >= 1){
                        contigID[ii] = contigID[ii]
                       #chrom = unique(dat1$RefContigID)
                        #chrom2 = unique(dat1$RefContigID)
                        
                        #modmolCmapdat <- modmolcmap(molcmap = moldat, xmapdata = xmapdat, cntNick = cntNick, cntMeth = cntMeth)
                        'if(length(chrom2) > 1){
                            nickrefLoc <- data.frame()
                            for(oo in 1:length(chrom2)){
                                nickrefLoc1 <- nickReference(refxmap = contigRefxmap, 
                                    refcmap = refcmap,
                                    contigcmap = querycmap,
                                    contigID = contigID[ii],
                                    chrom = chrom2[oo],
                                    returnMethod = "dataFrame")
                                nickrefLoc <- rbind(nickrefLoc, nickrefLoc1,row.names = NULL)
                            }
                        }else{
                            nickrefLoc <- nickReference(refxmap = contigRefxmap, 
                                refcmap = refcmap,
                                contigcmap = querycmap,
                                contigID = contigID[ii],
                                chrom = chrom2,
                                returnMethod = "dataFrame")
                        }'
						nickrefLoc <- nickReference(refxmap = contigRefxmap, 
                                refcmap = refcmap,
                                contigcmap = querycmap,
                                contigID = contigID[ii],
                                chrom = chrom[uu],
                                returnMethod = "dataFrame")
                        modmolCmapdat <- modmolcmap(molcmap = moldat, xmapdata = xmapdat, cntNick = cntNick, cntMeth = cntMeth)
                        methylNickLoc <- molNickLoc(nickrefLoc = nickrefLoc,
                            molcmap = modmolCmapdat,
                            xmapdata = xmapdat, 
                            contigID = contigID[ii],
							cntNick = cntNick, 
							cntMeth = cntMeth)
                        molcmapFinal <- molMethylLoc(methylNickLoc, cntNick = cntNick, cntMeth = cntMeth, outputType = c("dataframe"))
                        molcmapFinal1 <- molMethylRefLoc(molcmapFinal, cntNick = cntNick, cntMeth = cntMeth, outputType = c("dataframe"))
                        print (outpath)
						if(outputFormat == "individualText"){
                            write.table(molcmapFinal1, file.path(outpath, paste("modifiedMol_", contigID[ii], ".txt", sep = "")))
                        }else if(outputFormat == "dataFrame"){
                            modmol <- rbind(modmol,molcmapFinal1, row.names = NULL)
                        }else if(outputFormat == "multipleText"){
                            modmol <- rbind(modmol,molcmapFinal1, row.names = NULL)
                            print("Will print outside the loop")
                        }else{stop("outputMethod Incorrect!!!!")}
                    }else{ print(paste("Contig:", contigID[ii], " is absent."))
                        next;}
                    
                }
            }else{
			    chrom = chrom1
		        datcont <- contigRefxmap[which(contigRefxmap$RefContigID == chrom),]
			    contigID <- datcont$QryContigID
			    pa <- paste("---Processing Chromosome", 
			        paste("chr_", chrom, sep = ""), 
				    "--", sep = "")
			    print(pa)
			    log_print(pa)
			    pa1 <- paste("---Number of contigs for chromosme ", 
			        paste("chr_", chrom, sep = ""), ":", 
			    	length(contigID), "---", sep = "")
			
			    #pp1 <- paste("---Number of contigs for chromosme1:", length(contigID))
			    print(pa1)
			    log_print(pa1)
				pp6 <- (paste("Contig being analyzed:", contigID, sep = ""))
				print(pp6)
				log_print(pp6)
				if(dir.exists(file.path(path, paste("chr_", chrom, sep = "")))
			    & length(list.files(path = file.path(path, paste("chr_", chrom, sep = "")))) > 0){
			        pa7 <- paste("chr_", chrom,": alreay processed", sep = "")
				    log_print(pa7)
				next;
			    }else if(!dir.exists(file.path(path, paste("chr_", chrom, sep = "")))){
			        pp2 <- paste("---Creating directory for Chromosome", chrom, "--", sep = "")
			        print(paste("---Creating directory for Chromosome", chrom, "--", sep = ""))
			        log_print(pp2)
			        dir.create(file.path(path, paste("chr_", chrom, sep = "")))
				    outpath <- file.path(path, paste("chr_", chrom, sep = ""))
			    }else{
			        outpath <- file.path(path, paste("chr_", chrom, sep = ""))
			    }
                'dat1 <- refxmap[which(refxmap$QryContigID == contigID),]
                chrom = unique(dat1$RefContigID)'
                contigmap <- list.files (cmaploc_mol, pattern = paste("exp_refineFinal1_contig",contigID,"_r.cmap",sep = ""), full.names = TRUE)
                molmap <- list.files (cmaploc_mol, pattern = paste("exp_refineFinal1_contig",contigID,"_q.cmap",sep = ""), full.names = TRUE)
                xmap <- list.files (cmaploc_mol, pattern = paste("exp_refineFinal1_contig",contigID,".xmap",sep = ""), full.names = TRUE)
                #contigdat <- readingCmap(contigmap)
                if(molmap == NA){
                    print(paste("exp_refineFinal1_contig",contigID,"_q.cmap", "is discarded from analysis!!!"))
                    pp4 <- paste("exp_refineFinal1_contig",contigID,"_q.cmap", "is discarded from analysis!!!")
					log_print(pp4)
                    next;
                }else{
                    moldat <- readingCmap(molmap)
                }
                xmapdat <- readingXmap(xmap = xmap)
		    	if(is.na(moldat) == TRUE 
		    		    | is.na(xmapdat) == TRUE){
						log_print("Either xmap,molecule cmap or reference cmap file empty !!!!")
		    			print("Either xmap,molecule cmap or reference cmap file empty !!!!")
		    			next;
		    		}else{moldat = moldat; xmapdat = xmapdat}
                #dat1 <- contigRefxmap[which(contigRefxmap$QryContigID == contigID[ii]),]
                if(nrow(dat1) >= 1){
                        contigID = contigID
                       #chrom = unique(dat1$RefContigID)
                        chrom2 = unique(dat1$RefContigID)
                        xmapdat <- readingXmap(xmap)
                        #modmolCmapdat <- modmolcmap(molcmap = moldat, xmapdata = xmapdat, cntNick = cntNick, cntMeth = cntMeth)
                        if(length(chrom2) > 1){
                            nickrefLoc <- data.frame()
                            for(oo in 1:length(chrom2)){
                                nickrefLoc1 <- nickReference(refxmap = contigRefxmap, 
                                    refcmap = refcmap,
                                    contigcmap = querycmap,
                                    contigID = contigID,
                                    chrom = chrom2[oo],
                                    returnMethod = "dataFrame")
                                nickrefLoc <- rbind(nickrefLoc, nickrefLoc1,row.names = NULL)
                            }
                        }else{
                            nickrefLoc <- nickReference(refxmap = contigRefxmap, 
                                refcmap = refcmap,
                                contigcmap = querycmap,
                                contigID = contigID,
                                chrom = chrom2,
                                returnMethod = "dataFrame")
                        }
                        modmolCmapdat <- modmolcmap(molcmap = moldat, xmapdata = xmapdat, cntNick = cntNick, cntMeth = cntMeth)
                        methylNickLoc <- molNickLoc(nickrefLoc = nickrefLoc,
                            molcmap = modmolCmapdat,
                            xmapdata = xmapdat, 
                            contigID = contigID,
							cntNick = cntNick, 
							cntMeth = cntMeth)
                        molcmapFinal <- molMethylLoc(methylNickLoc, cntNick = cntNick, cntMeth = cntMeth, outputType = c("dataframe"))
                        molcmapFinal1 <- molMethylRefLoc(molcmapFinal, cntNick = cntNick, cntMeth = cntMeth, outputType = c("dataframe"))
						print (outpath)
                        if(outputFormat == "individualText"){
                            write.table(molcmapFinal1, file.path(outpath, paste("modifiedMol_", contigID, ".txt", sep = "")))
                        }else if(outputFormat == "dataFrame"){
                            modmol <- rbind(modmol,molcmapFinal1,row.names = NULL)
                        }else if(outputFormat == "multipleText"){
                            modmol <- rbind(modmol,molcmapFinal1,row.names = NULL)
                            print("Will print outside the loop")
                        }else{stop("outputMethod Incorrect!!!!")}
                    }else{ print(paste("Contig:", contigID, " is absent."))
                        next;}
                    
        }
        if(outputFormat == "multipleText"){
            if(length(unique(modmol$ContigID)) > 1){
                contig <- unique(modmol$ContigID)
                contigname <- paste(contig, collapse = "_")
                write.table(modmol, file.path(outpath, paste("modifiedMol_", contigname, ".txt", sep = "")))
            }else if (length(unique(modmol$ContigID)) == 1){
            print("Individual contig files have been written !!")
            }else{stop("No contigs found !!!")}
        }else if(outputFormat == "individualText"){
            log_print("Individual contig files have already been written !!")
        }else if(outputFormat == "dataFrame"){
            return(modmol)
        }else{stop("outputMethod Incorrect!!!!")}
	    pp6 <- paste("---Finished Processing", chrom[uu], "--", sep = "")
		log_print(pp6)
	}
        
    }else if(selectedRegion == TRUE){
	    print("Selected region is TRUE")
		dir.create(file.path(path, "Selected"))
        chrom = chrom;
        startPos = startPos; endPos = endPos;
        dat1 <- contigextract(contigRefxmap = contigRefxmap,
            chrom = chrom, startPos = startPos, endPos = endPos)
        if(nrow(dat1) >= 1){
            contigID <- dat1$contigID 
            chrom <- unique(dat1$chrom)
        }else{stop("Selected region does not have a contig mapped to it!!!")}
        print(paste("Total number of contigs to be processed:", length(contigID)))
		'pp7 <- paste("---Creating directory for Chromosome", chrom, "--", sep = "")
		print(paste("---Creating directory for Chromosome", chrom, "--", sep = ""))
		outpath <- dir.create(file.path(path, paste("chr_", chrom[uu], sep = "")))'
		
        if(dir.exists(file.path(path, "Selected", paste("chr_", chrom, sep = "")))
	        & length(list.files(path = file.path(path, paste("chr_", "Selected", chrom, sep = "")))) > 0){
			    pa7 <- paste("chr_", chrom,": alreay processed", sep = "")
				log_print(pa7)
				
		}else if(!dir.exists(file.path(path, "Selected", paste("chr_", chrom, sep = "")))){
		    pp2 <- paste("---Creating directory for Chromosome", chrom, "--", sep = "")
		    print(paste("---Creating directory for Chromosome", chrom, "--", sep = ""))
		    log_print(pp2)
		    dir.create(file.path(path, "Selected", paste("chr_", chrom, sep = "")))
			outpath <- file.path(path, "Selected", paste("chr_", chrom, sep = ""))
		}else{
		    outpath <- file.path(path, "Selected", paste("chr_", chrom, sep = ""))
		}    
        if(length(contigID) > 1){
                print(nrow(contigRefxmap)); 
				print(nrow(refcmap));
                for (ii in 1:length(contigID)){
		    	#for (ii in 803:length(contigID)){
		    	    print(paste("Contig being analyzed:", contigID[ii], sep = ""))
					pp5 <- paste("Contig being analyzed:", contigID[ii], sep = "")
					log_print(pp5)
                    print("Reading contig, molecule and reference xmap")
					contigmap <- list.files (contigMolPath, pattern = paste("exp_refineFinal1_contig",contigID[ii],"_r.cmap",sep = ""), full.names = TRUE)
                    molmap <- list.files (contigMolPath, pattern = paste("exp_refineFinal1_contig",contigID[ii],"_q.cmap",sep = ""), full.names = TRUE)
                    xmap <- list.files (contigMolPath, pattern = paste("exp_refineFinal1_contig",contigID[ii],".xmap",sep = ""), full.names = TRUE)
                    #contigdat <- readingCmap(contigmap)
                    if(molmap == ""){
                        print(paste("exp_refineFinal1_contig",contigID[ii],"_q.cmap", "is discarded from analysis!!!"))
                        pp4 <- paste("exp_refineFinal1_contig",contigID[ii],"_q.cmap", "is discarded from analysis!!!")
						log_print(pp4)
						next;
                    }else{
                        moldat <- readingCmap(molmap)
                    }
                    dat1 <- contigRefxmap[which(contigRefxmap$QryContigID == contigID[ii]),]
                    xmapdat <- readingXmap(xmap)
                    if(all(is.na((moldat[1:3,]))) == TRUE 
		    		    | all(is.na((xmapdat[1:3,]))) == TRUE  
		    			| all(is.na((querycmap[1:3,]))) == TRUE){
						log_print("Either xmap,molecule cmap or reference cmap file empty !!!!")
		    			print("Either xmap,molecule cmap or reference cmap file empty !!!!")
		    			next;
		    		}else{moldat = moldat; xmapdat = xmapdat}
                    if(nrow(dat1) >= 1){
                        contigID[ii] = contigID[ii]
                       #chrom = unique(dat1$RefContigID)
                        chrom2 = unique(dat1$RefContigID)
                        
                        #modmolCmapdat <- modmolcmap(molcmap = moldat, xmapdata = xmapdat, cntNick = cntNick, cntMeth = cntMeth)
                        if(length(chrom2) > 1){
                            nickrefLoc <- data.frame()
                            for(oo in 1:length(chrom2)){
                                nickrefLoc1 <- nickReference(refxmap = contigRefxmap, 
                                    refcmap = refcmap,
                                    contigcmap = querycmap,
                                    contigID = contigID[ii],
                                    chrom = chrom2[oo],
                                    returnMethod = "dataFrame")
                                nickrefLoc <- rbind(nickrefLoc, nickrefLoc1,row.names = NULL)
                            }
                        }else{
                            nickrefLoc <- nickReference(refxmap = contigRefxmap, 
                                refcmap = refcmap,
                                contigcmap = querycmap,
                                contigID = contigID[ii],
                                chrom = chrom2,
                                returnMethod = "dataFrame")
                        }
                        modmolCmapdat <- modmolcmap(molcmap = moldat, xmapdata = xmapdat, cntNick = cntNick, cntMeth = cntMeth)
                        methylNickLoc <- molNickLoc(nickrefLoc = nickrefLoc,
                            molcmap = modmolCmapdat,
                            xmapdata = xmapdat, 
                            contigID = contigID[ii],
							cntNick = cntNick, 
							cntMeth = cntMeth)
                        molcmapFinal <- molMethylLoc(methylNickLoc, cntNick = cntNick, cntMeth = cntMeth, outputType = c("dataframe"))
                        molcmapFinal1 <- molMethylRefLoc(molcmapFinal, cntNick = cntNick, cntMeth = cntMeth, outputType = c("dataframe"))
                        if(outputFormat == "individualText"){
                            write.table(molcmapFinal1, file.path(outpath, paste("modifiedMol_", contigID[ii], ".txt", sep = "")))
                        }else if(outputFormat == "dataFrame"){
                            modmol <- rbind(modmol,molcmapFinal1, row.names = NULL)
                        }else if(outputFormat == "multipleText"){
                            modmol <- rbind(modmol,molcmapFinal1, row.names = NULL)
                            print("Will print outside the loop")
                        }else{stop("outputMethod Incorrect!!!!")}
                    }else{ print(paste("Contig:", contigID[ii], " is absent."))
                        next;}
                    
                }
            }else{
		        print(paste("Contig being analyzed:", contigID, sep = ""))
                'dat1 <- refxmap[which(refxmap$QryContigID == contigID),]
                chrom = unique(dat1$RefContigID)'
                contigmap <- list.files (cmaploc_mol, pattern = paste("exp_refineFinal1_contig",contigID,"_r.cmap",sep = ""), full.names = TRUE)
                molmap <- list.files (cmaploc_mol, pattern = paste("exp_refineFinal1_contig",contigID,"_q.cmap",sep = ""), full.names = TRUE)
                xmap <- list.files (cmaploc_mol, pattern = paste("exp_refineFinal1_contig",contigID,".xmap",sep = ""), full.names = TRUE)
                #contigdat <- readingCmap(contigmap)
                if(molmap == NA){
                    print(paste("exp_refineFinal1_contig",contigID,"_q.cmap", "is discarded from analysis!!!"))
                    pp4 <- paste("exp_refineFinal1_contig",contigID,"_q.cmap", "is discarded from analysis!!!")
					log_print(pp4)
                    
                }else{
                    moldat <- readingCmap(molmap)
                }
                xmapdat <- readingXmap(xmap = xmap)
		    	if(is.na(moldat) == TRUE 
		    		    | is.na(xmapdat) == TRUE){
						log_print("Either xmap,molecule cmap or reference cmap file empty !!!!")
		    			print("Either xmap,molecule cmap or reference cmap file empty !!!!")
		    			
		    		}else{moldat = moldat; xmapdat = xmapdat}
                #dat1 <- contigRefxmap[which(contigRefxmap$QryContigID == contigID[ii]),]
                if(nrow(dat1) >= 1){
                        contigID = contigID
                       #chrom = unique(dat1$RefContigID)
                        chrom2 = unique(dat1$RefContigID)
                        xmapdat <- readingXmap(xmap)
                        #modmolCmapdat <- modmolcmap(molcmap = moldat, xmapdata = xmapdat, cntNick = cntNick, cntMeth = cntMeth)
                        if(length(chrom2) > 1){
                            nickrefLoc <- data.frame()
                            for(oo in 1:length(chrom2)){
                                nickrefLoc1 <- nickReference(refxmap = contigRefxmap, 
                                    refcmap = refcmap,
                                    contigcmap = querycmap,
                                    contigID = contigID,
                                    chrom = chrom2[oo],
                                    returnMethod = "dataFrame")
                                nickrefLoc <- rbind(nickrefLoc, nickrefLoc1,row.names = NULL)
                            }
                        }else{
                            nickrefLoc <- nickReference(refxmap = contigRefxmap, 
                                refcmap = refcmap,
                                contigcmap = querycmap,
                                contigID = contigID,
                                chrom = chrom2,
                                returnMethod = "dataFrame")
                        }
                        modmolCmapdat <- modmolcmap(molcmap = moldat, xmapdata = xmapdat, cntNick = cntNick, cntMeth = cntMeth)
                        methylNickLoc <- molNickLoc(nickrefLoc = nickrefLoc,
                            molcmap = modmolCmapdat,
                            xmapdata = xmapdat, 
                            contigID = contigID,
							cntNick = cntNick, 
							cntMeth = cntMeth)
                        molcmapFinal <- molMethylLoc(methylNickLoc, cntNick = cntNick, cntMeth = cntMeth, outputType = c("dataframe"))
                        molcmapFinal1 <- molMethylRefLoc(molcmapFinal, cntNick = cntNick, cntMeth = cntMeth, outputType = c("dataframe"))
                        if(outputFormat == "individualText"){
                            write.table(molcmapFinal1, file.path(outpath, paste("modifiedMol_", contigID, ".txt", sep = "")))
                        }else if(outputFormat == "dataFrame"){
                            modmol <- rbind(modmol,molcmapFinal1,row.names = NULL)
                        }else if(outputFormat == "multipleText"){
                            modmol <- rbind(modmol,molcmapFinal1,row.names = NULL)
                            print("Will print outside the loop")
                        }else{stop("outputMethod Incorrect!!!!")}
                    }else{ print(paste("Contig:", contigID, " is absent."))
                        }
                    
        }
        if(outputFormat == "multipleText"){
            if(length(unique(modmol$ContigID)) > 1){
                contig <- unique(modmol$ContigID)
                contigname <- paste(contig, collapse = "_")
                write.table(modmol, file.path(outpath, paste("modifiedMol_", contigname, ".txt", sep = "")))
            }else if (length(unique(modmol$ContigID)) == 1){
            print("Individual contig files have been written !!")
            }else{stop("No contigs found !!!")}
        }else if(outputFormat == "individualText"){
            log_print("Individual contig files have already been written !!")
        }else if(outputFormat == "dataFrame"){
            return(modmol)
        }else{stop("outputMethod Incorrect!!!!")}
        pp10 <- paste("---Finished Processing Selected Region---")
		log_print(pp10)
    
    
    }else{stop("Selected Region Not Selected!!!")}
    
    log_print("Molecule cmaps modified !!!")
    print("Molecule cmaps modified !!!")
	log_print("---Closing Log Files!!!---")
	log_close()
    #return(paste("query methylation Labels Modified"))
    #return (fname1)
}




