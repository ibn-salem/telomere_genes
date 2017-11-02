

#' Add distance in bp to next telomere
#' 
add_telomere_dist <- function(gr){
  
  stopifnot(length(seqinfo(gr)) > 1)
  
  # build GRanges object for human telomeres
  chromGR <- GRanges(seqinfo(gr))
  chromStartsGR <- flank(chromGR, width = -1, start = TRUE)
  chromStartsGR$telomere <- paste0(seqnames(chromStartsGR), "_start")
  
  chromEndsGR <- flank(chromGR, width = -1, start = FALSE)
  chromEndsGR$telomere <- paste0(seqnames(chromEndsGR), "_end")
  
  telomereGR <- c(chromStartsGR, chromEndsGR)
  
  # get dist to next telomere
  hits <- distanceToNearest(gr, telomereGR)
  
  # assume all genes accoure only once
  stopifnot(identical(queryHits(hits), seq(length(gr))))
  
  # build data.frame with telomere and distance
  telomereDF <- hits %>% 
    as.data.frame() %>% 
    as.tibble() %>% 
    mutate(
      telomere = telomereGR$telomere[subjectHits],
      telomere_distance = distance
    ) %>% 
    select(telomere, telomere_distance)
  
  # add telomere annotation to gr
  gr$telomere = telomereDF$telomere
  gr$telomere_dist = telomereDF$telomere_distance
  
  return(gr)
}
