get450KProbeMapping <- function(probeIDs, platform='HM450', genome='hg19'){
  hm450 <- FDb.InfiniumMethylation.hg19::getPlatform(platform = platform, genome = genome)
  probes <- hm450[probeIDs]

  TSS = FDb.InfiniumMethylation.hg19::getNearestTSS(probes)
  TSS = rownameToFirstColumn(TSS,'methProbeIDs')
  TSS = dplyr::select(TSS, methProbeIDs, distance, nearestGeneSymbol, nearestTranscript)
  setnames(TSS,c("distance", "nearestGeneSymbol", "nearestTranscript"),
           c("distanceToTSS", "nearestTSS", "nearestTSS.ID"))

  Tx = FDb.InfiniumMethylation.hg19::getNearestTranscript(probes)
  Tx = rownameToFirstColumn(Tx, 'methProbeIDs')
  Tx = dplyr::select(Tx, methProbeIDs, distance, nearestGeneSymbol, nearestTranscript)
  setnames(Tx,c("distance", "nearestGeneSymbol", "nearestTranscript"),
           c("distanceToTx", "nearestTx", "nearestTx.ID"))

  hm450 = rownameToFirstColumn(hm450, 'methProbeIDs')

  Annotation = plyr::join_all(list(hm450,TSS,Tx), by = 'methProbeIDs', match = 'all')

  return(list(Annotation = Annotation))
}
