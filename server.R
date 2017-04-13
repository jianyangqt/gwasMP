######################################
# GWAS Mapping Precision (gwasMP)
# @Author zhili zheng <zhilizheng@uq.edu.au>;
# @Author Yang Wu <y.wu2@uq.edu.au>;
# @Author Jian Yang <jian.yang@uq.edu.au>;
# @Bugs report to Jian Yang <jian.yang@uq.edu.au>
# @Updated Feb.26 2017
# @Licence GPL v3
######################################

source("helpers.R")
mafbin=c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)

load("data.rdata")
data = common

dthresh=1e7
pthresh=5e-8
pdthreshold=seq(5,100,by=5)
ldthreshold=seq(0,1,by=0.05)
data_type = c("wgs","1kg3","1kg1","hap")
legends = c("WGS","1KGP3","1KGP1","HAPMAP2")

shinyServer(
  function(input, output) {
      dataInput <- reactive({    
            data = switch(input$maf,
                "Common (MAF > 0.01)" = common,
				"0.01 < MAF < 0.05" = common[(common$MAF_qtl>=0.01) & (common$MAF_qtl<0.05), ],
				"0.05 < MAF < 0.10" = common[(common$MAF_qtl>=0.05) & (common$MAF_qtl<0.10), ],
				"0.10 < MAF < 0.20" = common[(common$MAF_qtl>=0.10) & (common$MAF_qtl<0.20), ],
				"0.20 < MAF < 0.30" = common[(common$MAF_qtl>=0.20) & (common$MAF_qtl<0.30), ],
				"0.30 < MAF < 0.40" = common[(common$MAF_qtl>=0.30) & (common$MAF_qtl<0.40), ],
				"0.40 < MAF < 0.50" = common[(common$MAF_qtl>=0.40) & (common$MAF_qtl<0.50), ],
                "Rare (0.0004 < MAF < 0.01)" = rare[(rare$MAF_qtl>=0.0004) & (rare$MAF_qtl<0.01),] 
			)
				
			data_list = list()
			data_list[["wgs"]]=which(data$pval_uk10k<=pthresh & data$distance_wgs<=dthresh)
            data_list[["1kg3"]]=which(data$pval_1kg3<=pthresh & data$distance_1kg3<=dthresh)
            data_list[["1kg1"]]=which(data$pval_1kg1<=pthresh & data$distance_1kg1<=dthresh)
			if(!all(is.na(data[1,"MAF_hap"]))){
				data_list[["hap"]]=which(data$pval_hap<=pthresh & data$distance_hap<=dthresh)
				data_list[["has_hap"]]=TRUE
				mafthreshold = seq(0,0.05,by=0.0025)
			}else{
				data_list[["hap"]]=c()
				data_list[["has_hap"]]=FALSE
				mafthreshold = seq(0,0.003,by=0.00015)
			}
            
			data_list[["data"]] = data

			
			pd_maps = data.frame(mark=c(),dis=c(),prop=c(),stringsAsFactors=FALSE)
			LD_maps = data.frame(mark=c(),ld=c(),prop=c(),stringsAsFactors=FALSE)
			maf_maps = data.frame(mark=c(),maf=c(),prop=c(),stringsAsFactors=FALSE)
			
			for(mark in data_type){
				data_index = data_list[[mark]]
				cur_data = data_list[["data"]][data_index,]
				data_dis = cur_data[,paste0("distance_",mark)] / 1000
				data_LD = cur_data[, paste0("rsq_",mark)]
				data_maf = abs(cur_data[, paste0("MAF_",mark)] - cur_data$MAF_qtl)
				pd_maps = rbind(pd_maps, data.frame(mark=mark,dis=pdthreshold,prop=propfunc(pdthreshold,data_dis,0)))
				LD_maps = rbind(LD_maps, data.frame(mark=mark,ld=ldthreshold,prop=propfunc(ldthreshold,data_LD,1)))
				maf_maps = rbind(maf_maps, data.frame(mark=mark,maf=mafthreshold,prop=propfunc(mafthreshold,data_maf,0)))
			}
			
			data_list[["pd_maps"]] = pd_maps
			data_list[["LD_maps"]] = LD_maps
			data_list[["maf_maps"]] = maf_maps
			return(data_list)
			
      })

      cacPropData <- reactive({
			validate(
				need(input$ld>=0 & input$ld<=1.0, 'LD threshold ranges from 0 to 1'),
				need(input$pd>=0 & input$pd<=1000, 'Physical distance ranges from 0 to 1000 Kb'),
				need(input$maf_diff>=0 & input$maf_diff<=0.05, 'MAF difference ranges from 0 to 0.05')
			)
			data_list = dataInput()
            LD = input$ld
			pd = input$pd
			maf = input$maf_diff

			re_props = data.frame(mark=c(),type=c(),value=c(),prop=c(),stringsAsFactors=FALSE)
			for(mark in data_type){
				props = match_ld_pos(LD, pd, maf, data_list, mark)
				re_props = rbind(re_props, data.frame(mark=mark, type="LD", value=as.character(round(LD,2)), prop=props[["LD_prop"]]))
				re_props = rbind(re_props, data.frame(mark=mark, type="pos", value=as.character(round(pd,2)), prop=props[["pos_prop"]]))
				re_props = rbind(re_props, data.frame(mark=mark, type="maf", value=as.character(round(maf,4)), prop=props[["maf_prop"]]))
				re_props = rbind(re_props, data.frame(mark=mark, type="both", value=paste("LD r<sup>2</sup> > ",LD, " & Distance < ",pd," Kb", " & MAF_diff < ",maf,sep=""), prop=props[["both_prop"]]))
			}
			
			list(re_props=re_props, LD=LD, pd=pd, maf = maf,data=data_list)
      })


      output$plot <- renderPlot({
            ############################################ plot figures ##############################################
            data_list = cacPropData()
			
            par(mar=c(4,5,1,2), mfrow=c(1,3))
			           		
			maps = data_list[["data"]][["pd_maps"]]
			pd = data_list[["pd"]]
			pd_max = max(pdthreshold)
			if(pd > pd_max){
				new_pdthresholds = seq(pd_max+10, ceiling((pd-pd_max)/20)*20+100, 10)
				for(mark in data_type){
					data_index = data_list[["data"]][[mark]]
					cur_data = data_list[["data"]][["data"]][data_index,]
					data_dis = cur_data[,paste0("distance_",mark)] / 1000
					maps = rbind(maps, data.frame(mark=mark,dis=new_pdthresholds,prop=propfunc(new_pdthresholds,data_dis,0)))
				}
				max_dis = max(new_pdthresholds)
			}else{
				max_dis = max(pdthreshold)
			}
            plot(maps[maps$mark=="wgs","dis"],maps[maps$mark=="wgs","prop"],type="b",col="orange",pch=23,lty=1,cex=0.8,
			    cex.axis=1.6,cex.lab=1.6,lwd=1.5,xlab="Distance (Kb)",ylab="Proportion (%)",bty="n",font.main=1,ylim=c(0,100),xlim=c(0,max_dis))
            lines(maps[maps$mark=="1kg3","dis"],maps[maps$mark=="1kg3","prop"],type="b",col="blue",pch=17,lty=1,cex=0.8,lwd=2)
            lines(maps[maps$mark=="1kg1","dis"],maps[maps$mark=="1kg1","prop"],type="b",col="black",pch=16,lty=1,cex=0.8,lwd=2)
            lines(maps[maps$mark=="hap","dis"],maps[maps$mark=="hap","prop"],type="b",col="green",pch=24,lty=1,cex=0.8,lwd=2)
            abline(v=pd,lty=2,col="red")
			if(data_list[["data"]][["has_hap"]]){
				 legend("bottomright",legends,lty=c(1,1,1,1),pch=c(23,17,16,24),col=c("orange","blue","black","green"),bty="n",cex=1.3,pt.cex=1)
			}else{
				 legend("bottomright",legends[1:3],lty=c(1,1,1),pch=c(23,17,16),col=c("orange","blue","black"),bty="n",cex=1.3,pt.cex=1)
			}

			maps = data_list[["data"]][["LD_maps"]]
			LD = data_list[["LD"]]
			plot(maps[maps$mark=="wgs","ld"],maps[maps$mark=="wgs","prop"],type="b",col="orange",pch=23,lty=1,cex=0.8,
			    cex.axis=1.6,cex.lab=1.6,lwd=2,xlab=expression(paste("LD (",r^2,")",sep="")), ylab="Proportion (%)",bty="n",ylim=c(0,100),xlim=c(0,1))
            lines(maps[maps$mark=="1kg3","ld"],maps[maps$mark=="1kg3","prop"],type="b",col="blue",pch=17,lty=1,cex=0.8,lwd=2)
            lines(maps[maps$mark=="1kg1","ld"],maps[maps$mark=="1kg1","prop"],type="b",col="black",pch=16,lty=1,cex=0.8,lwd=2)
            lines(maps[maps$mark=="hap","ld"],maps[maps$mark=="hap","prop"],type="b",col="green",pch=24,lty=1,cex=0.8,lwd=2)
            abline(v=LD,lty=2,col="red")
			if(data_list[["data"]][["has_hap"]]){
			    legend("bottomleft",legends,lty=c(1,1,1,1),pch=c(23,17,16,24),col=c("orange","blue","black","green"),bty="n",cex=1.3,pt.cex=1)
			}else{
				 legend("bottomleft",legends[1:3],lty=c(1,1,1),pch=c(23,17,16),col=c("orange","blue","black"),bty="n",cex=1.3,pt.cex=1)
			}

			maps = data_list[["data"]][["maf_maps"]]
			pd = data_list[["maf"]]
            plot(maps[maps$mark=="wgs","maf"],maps[maps$mark=="wgs","prop"],type="b",col="orange",pch=23,lty=1,cex=0.8,cex.axis=1.6,
			    cex.lab=1.6,lwd=2,xlab="MAF difference",ylab="Proportion (%)",bty="n",ylim=c(0,100),xlim=range(maps[maps$mark=="wgs","maf"]))
            lines(maps[maps$mark=="1kg3","maf"],maps[maps$mark=="1kg3","prop"],type="b",col="blue",pch=17,lty=1,cex=0.8,lwd=2)
            lines(maps[maps$mark=="1kg1","maf"],maps[maps$mark=="1kg1","prop"],type="b",col="black",pch=16,lty=1,cex=0.8,lwd=2)
            lines(maps[maps$mark=="hap","maf"],maps[maps$mark=="hap","prop"],type="b",col="green",pch=24,lty=1,cex=0.8,lwd=2)
            abline(v=pd,lty=2,col="red")
			if(data_list[["data"]][["has_hap"]]){
				 legend("bottomright",legends,lty=c(1,1,1,1),pch=c(23,17,16,24),col=c("orange","blue","black","green"),bty="n",cex=1.3,pt.cex=1)
			}else{
				 legend("bottomright",legends[1:3],lty=c(1,1,1),pch=c(23,17,16),col=c("orange","blue","black"),bty="n",cex=1.3,pt.cex=1)
			}

			
      })

      output$prop <- renderTable({
		re_props = (cacPropData())[["re_props"]]
		pd_prop = as.matrix(re_props[re_props$type=="pos",c(3,4)])
		LD_prop = as.matrix(re_props[re_props$type=="LD",c(3,4)])
		maf_prop = as.matrix(re_props[re_props$type=="maf",c(3,4)])
		sep = matrix(rep(" ",4),nrow=4)
		prop = cbind(pd_prop, sep, LD_prop, sep, maf_prop)
		colnames(prop) = c("Distance (Kb)", "Proportion (%)", "|", "LD (r<sup>2</sup>)", "Proportion (%)", "|", "MAF_diff", "Proportion (%)")
		rownames(prop) = legends
		prop
	  }, rownames=TRUE, align="lcccccccc", sanitize.text.function = function(x) x)
	  

	  output$both_prop <- renderTable({
		 re_props = (cacPropData())[["re_props"]]
		 both_prop = as.matrix(re_props[re_props$type=="both",c(3,4)])
		 colnames(both_prop) = c("All thresholds", "Proportion (%)")
		 rownames(both_prop) = legends
		 both_prop
	  }, rownames=TRUE, align="lcc", sanitize.text.function = function(x) x)
  }
)