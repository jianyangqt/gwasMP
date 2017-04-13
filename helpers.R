######################################
# GWAS Mapping Precision (gwasMP)
# @Author zhili zheng <zhilizheng@uq.edu.au>;
# @Author Yang Wu <y.wu2@uq.edu.au>;
# @Author Jian Yang <jian.yang@uq.edu.au>;
# @Bugs report to Jian Yang <jian.yang@uq.edu.au>
# @Updated Feb.26 2017
# @Licence GPL v3
######################################

propfunc<-function(threshold,dataset,flag){
	prop=c()
	if(flag==0){
		for(i in 1:length(threshold)){
		prop[i]=(length(which(dataset<=threshold[i]))/length(dataset))*100
		}
	}
	if(flag==1){
		for(i in 1:length(threshold)){
		prop[i]=(length(which(dataset>=threshold[i]))/length(dataset))*100
		}
	}
	return(prop)
}

targetfunc<-function(prop,dataset){
	target=sort(dataset)[floor(length(dataset)*prop)];
	return(target)
}


match_ld_pos <- function(LD, pos, maf, data_list, mark){
	data_index = data_list[[mark]]
	data_set = data_list[["data"]][data_index,]

	LD_mark = data_set[[paste0("rsq_",mark)]] >= LD
	LD_prop = round(100 * sum(LD_mark, na.rm = TRUE) / dim(data_set)[1], 2)

	pos_mark = data_set[[paste0("distance_",mark)]] <= (pos * 1000)
	pos_prop = round(100 * sum(pos_mark, na.rm = TRUE) / dim(data_set)[1], 2)

	maf_mark = (abs(data_set[[paste0("MAF_",mark)]] - data_set$MAF_qtl) <= maf)
	maf_prop = round(100 * sum(maf_mark, na.rm = TRUE) / dim(data_set)[1], 2)

	both_mark = LD_mark & pos_mark & maf_mark
	both_prop = round(100 * sum(both_mark, na.rm = TRUE) /dim(data_set)[1], 2)

	return(list(LD_prop=LD_prop, pos_prop=pos_prop, maf_prop=maf_prop, both_prop=both_prop))

}