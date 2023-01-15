# Mann Whitney U tests

# Text comparisons
headings = c('name','t.test.p','mann.whitney.p','shapiro.wilks.p.x','shapiro.wilks.p.y','t.test.p.removed')
text_comparisons = data.frame(matrix(nrow = 0,ncol = length(headings)))
names(text_comparisons) = headings

# Variant modified vs 1
i=1
#t.test p-value
text_comparisons[i,]=write_stats('vsv on variants vs 1', log10(vsv_improvement_variant))
i=i+1
text_comparisons[i,]=write_stats('vsv on variants vs ancestral', log10(vsv_improvement_variant),log10(vsv_improvement_ancestral))
i=i+1
text_comparisons[i,]=write_stats('vsv matched vs nonmatched', log10(vsv_improvement_variant_matched),log10(vsv_improvement_variant_nonmatched))
i=i+1
text_comparisons[i,]=write_stats('monovalent vs bivalent',log10(vsv_rises_monovalent), log10(vsv_rises_bivalent))
i=i+1
text_comparisons[i,]=write_stats('uninf vs inf',log10(vsv_rises_prioruninf), log10(vsv_rises_priorinfect))
i=i+1
text_comparisons[i,]=write_stats('2 vs 3 doses',log10(vsv_rises_priorprimary), log10(vsv_rises_priorboosted))
i=i+1
text_comparisons[i,]=write_stats('severe above 0',effs_avged[which(start_effs==.5),which(col.names=='severe.improved.diff'),])




  
write_stats = function(name, x, y=NULL) {
  if (!is.null(y)){
    cd.y = cooks.distance(lm(y~1))
    ynew = y[cd.y<(4/length(x))]
  }
    
  cd.x = cooks.distance(lm(x~1))
  xnew = x[cd.x<(4/length(x))]
  # can only do shapiro test on  5000 elements 
  if (length(x)>5000){
    xshap = x[sample.int(length(x), 5000)]
  } else {xshap = x}
  if (length(y)>5000){
    yshap = y[sample.int(length(y), 5000)]
  } else {yshap = y}
  if (is.null(y)){
    output = c(name,t.test(x)$p.value, wilcox.test(x)$p.value, shapiro.test(xshap)$p.value, NA,t.test(xnew)$p.value)
  } else {
    output = c(name,t.test(x,y)$p.value, wilcox.test(x,y)$p.value, shapiro.test(xshap)$p.value, shapiro.test(yshap)$p.value,t.test(xnew,ynew)$p.value)
  }
  output
}