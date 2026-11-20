# Complete75 / phosphosite overlap and biological interpretation

The student explicitly requested the overlap investigation previously deferred.
Matched their12 proteins/13 distinct sites/14 site-medium hits by exact reference
UniProt accession and independently verified JE2 locus/protein_key. Three proteins
overlap75:IsaA(Q2FV52),Pyk(Q2FXM9),QoxB(Q2FZK0). Previous49hadIsaA/Pyk;67hadQoxB.
No enrichment test on these selected lists. All14hit selections reproduce the
step9 observed,FDR<.1,abslog2effect>1 exclusive thresholded sets. Exclusive does
not mean present only in one mutant or a significant difference between mutants.

NewR scripts 06_compare_phosphosite_overlap.R and
06_write_phosphosite_overlap_report.R save membership,source effects/FDR,
parent comparisons,JE2 residue checks,localization summaries,network evidence,
primary literature manifest and Markdown explanation. No raw per-sample table
was copied into reports; a rendered plot shows observed normalized samplevalues
for the3totalproteins/3sites,TSB/SNonly. Allsource data remain external/ignored.

TSBphosphosite results:IsaAS194stk/anc+3.04020,q.015658,parent+2.27354,q6.62e-10;
PykT537stp/anc+1.88484,q.006037,parent+.14244,q.014442;
QoxBT552stp/anc+1.69406,q.069911,parent+.07457,q.53014.
Allpositive mean higher signal/abundance in the indicated mutant than ancestor.
All3parents have3/3detections in relevantgroups,allinteractionmodelsobserved.
Noneof9correspondingphosphositeTSBvsSNinteractions passesFDR.1(minq.615981).
Do not call the sites PASN-dependent phosphorylation responses or direct
Stk1/Stp1substrates. No protein-adjusted occupancy test or newFDRcomputed.

Pyk18/18sitevalueslocalized.9279-.9574;IsaA16/18,observedLP0-.8964;
QoxB17/18,.4798-.9496. BestFragPipeprobabilities.9574/.8964/.9496,respectively.
IsaA0LP occurs with a finitequantifiedvalue andpositiveFragPipeintensity;
it is not just a missingquantification artifact. Preserve this uncertainty.
All3positions checkedagainstJE2FASTA. PykT537windowALGLITEENGITSPSAIVGLEK
matches corresponding peptide in Burnside2010Table6,WTNewmanandDelta-stk1.

Literature synthesis:IsaAcellwallhydrolase,contributes tobiofilm/betalactam
phenotypes;Stapleton2007 also links SrrAB toisaAexpression. Ma2022findspurgene/
proteinresponses toisaAdisruption,so mainfunctionalone neverexcluded indirect
purineassociation. Pykenergy/pyruvatesupplymakesrespiratorytoxinadaptation
plausible;Vasu2015PknBinvitrophosphorylation/activityworkdoesnotidentifyT537;
Burnside2010supportspossibilityofStk1-independentphosphorylation. QoxBterminal
menaquinoloxidasefunctionHammer2016andHQNOrespiratoryinhibitionRoman-Rodriguez
2025supportPASNcontext,butnotT552functionorHQNOasthesolecauseinourmixture.
Schwyter2025andPrust2021coursepapersprovidecontext,notdirectsitevalidation.

Plotvisuallyinspected:all6panels/18groupsperassaypresent,missingcountslabelled,
separateyscales,noSNpulse,noimputedpoints. Independentread-onlyagentverified
rawderivedeffects/FDRs,all14hits,none9siteinteractionssignificant,JE2positions,
andFragPipelocalizationfacts. Existingpackagebuild/locale warnings only.
No existing75networkfigure,ORA/GSEA,sourceDEAorhistoricalreportsmodified.

AIassistance:rootCodexwroteRanalysis/reportandcheckedQoxBliterature;twoagents
researchedIsaA/Pykprimarypapers;afourthindependentlyauditedexactoverlapand
sourcequantitative/localizationevidence. No commit/pushrequested.
