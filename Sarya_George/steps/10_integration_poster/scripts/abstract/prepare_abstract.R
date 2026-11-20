# Prepare the lecture-recipe draft, check word counts and retrieve bibliographic
# metadata from PubMed. No experimental estimates are changed or reanalysed.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/10_integration_poster/reports/abstract")
body <- paste(readLines(file.path(out,"abstract_body.txt"),warn=FALSE),collapse=" ")
title <- paste(readLines(file.path(out,"title.txt"),warn=FALSE),collapse=" ")
count_words <- function(x) length(strsplit(trimws(x),"[[:space:]]+")[[1]])
body_count <- count_words(body); title_count <- count_words(title)
stopifnot(body_count<=300L,body_count+title_count<=300L,
  !grepl("\n",body),grepl("Here we show",body,fixed=TRUE))
italic_species <- function(x) {
  for(species in c("Pseudomonas aeruginosa","Staphylococcus aureus")) {
    x <- gsub(species,paste0("*",species,"*"),x,fixed=TRUE)
  }
  x
}
writeLines(c(paste0("# ",italic_species(title)),"",italic_species(body),"",
  paste0("Abstract body: ",body_count," words; body plus title: ",body_count+title_count,
    " words (whitespace-delimited; hyphenated terms counted as one)."),"",
  "Working abstract title. Per the student's instruction, poster-title coordination is outside this draft. Author and affiliation details are not inferred.","",
  "This is a draft for student review. The body follows the lecture sequence: general background, specific background, knowledge gap, Here we show, results, implications. References and recipe/source notes are separate from the submission text."),
  file.path(out,"abstract_draft.md"))
write_tsv(tibble(component=c("abstract_body","title","body_plus_title"),
  words=c(body_count,title_count,body_count+title_count),limit=300L),file.path(out,"word_count.tsv"))

# Primary records verified in PubMed-domain literature/author searches. Preserve
# exact source metadata so Zotero/EndNote imports do not rely on invented authors.
selected <- tribble(
  ~pmid,~role,~expected_doi,
  "40741853","Core: JE2 supernatant-induced protein responses","10.1128/msystems.00262-25",
  "33444734","Core: Stk1/Stp1 phosphoproteomics","10.1074/mcp.RA120.002232",
  "36595572","Core: phosphoregulation and stress adaptation","10.1126/scisignal.abj8194",
  "37646506","Supporting: secreted-compound inhibition and evolutionary context","10.1128/mbio.03153-22",
  "41114581","Supporting: respiratory inhibition and fermentation","10.1128/jb.00101-25")
r <- .api_request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
  params=list(db="pubmed",id=paste(selected$pmid,collapse=","),rettype="medline",retmode="text"),
  service="PubMed",refresh=FALSE)
raw <- file.path(REPO_ROOT,"data/raw/abstract_literature")
dir.create(raw,recursive=TRUE,showWarnings=FALSE)
writeLines(r$body,file.path(raw,"pubmed_selected.medline"),useBytes=TRUE)
lines <- strsplit(r$body,"\n",fixed=TRUE)[[1]]
starts <- which(grepl("^PMID-",lines)); ends <- c(starts[-1]-1L,length(lines))
read_record <- function(x) {
  fields <- list(); previous <- NULL
  for(line in x) {
    if(grepl("^[A-Z0-9]{2,4}[ ]*- ",line)) {
      tag <- trimws(sub("-.*$","",line))
      value <- sub("^[A-Z0-9]{2,4}[ ]*- ","",line)
      fields[[tag]] <- c(fields[[tag]],value); previous <- tag
    } else if(grepl("^      ",line) && !is.null(previous)) {
      n <- length(fields[[previous]])
      fields[[previous]][n] <- paste(fields[[previous]][n],trimws(line))
    }
  }
  fields
}
records <- lapply(seq_along(starts),function(i) read_record(lines[starts[i]:ends[i]]))
metadata <- bind_rows(lapply(records,function(x) {
  doi <- sub(" \\[doi\\]$","",x$AID[grepl("\\[doi\\]$",x$AID)])
  stopifnot(length(doi)==1L)
  tibble(pmid=x$PMID,title=str_squish(x$TI),year=substr(x$DP[1],1,4),journal=x$JT,
    volume=if(length(x$VI)) x$VI[1] else NA_character_,
    issue=if(length(x$IP)) x$IP[1] else NA_character_,
    pages=if(length(x$PG)) x$PG[1] else NA_character_,doi=doi,authors=paste(x$FAU,collapse="; "),
    publication_types=paste(x$PT,collapse="; "),
    journal_article_url=paste0("https://doi.org/",doi),
    pubmed_url=paste0("https://pubmed.ncbi.nlm.nih.gov/",x$PMID,"/"))
}))
metadata <- selected |> left_join(metadata,by="pmid",relationship="one-to-one")
metadata <- metadata |> mutate(reference_type="Primary experimental research",
  abstract_core=startsWith(role,"Core:"))
stopifnot(nrow(metadata)==5L,!anyNA(metadata$title),
  all(tolower(metadata$doi)==tolower(metadata$expected_doi)),
  metadata$year[metadata$pmid=="36595572"]=="2023",
  !any(grepl("(^|; )Review($|;)",metadata$publication_types)),sum(metadata$abstract_core)==3L)
# Primary status also checked against each article's original experimental
# methods/results. A Journal Article tag alone does not establish primary status.
write_tsv(metadata,file.path(out,"bibliography.tsv"))
make_ris <- function(refs) unlist(lapply(seq_len(nrow(refs)),function(i) {
  z <- refs[i,]
  authors <- strsplit(z$authors,"; ",fixed=TRUE)[[1]]
  c("TY  - JOUR",paste0("TI  - ",z$title),paste0("AU  - ",authors),
    paste0("PY  - ",z$year),paste0("JO  - ",z$journal),paste0("VL  - ",z$volume),
    if(!is.na(z$issue)) paste0("IS  - ",z$issue),
    if(!is.na(z$pages)) paste0("SP  - ",z$pages),
    paste0("DO  - ",z$doi),paste0("AN  - ",z$pmid),paste0("UR  - ",z$journal_article_url),"ER  - ","")
}),use.names=FALSE)
writeLines(make_ris(metadata),file.path(out,"references.ris"),useBytes=TRUE)
core <- filter(metadata,abstract_core)
writeLines(make_ris(core),file.path(out,"abstract_core_references.ris"),useBytes=TRUE)
citations <- vapply(seq_len(nrow(core)),function(i) {
  z <- core[i,]
  lead_author <- sub(",.*$","",strsplit(z$authors,"; ",fixed=TRUE)[[1]][1])
  paste0(i,". ",lead_author," et al. (",z$year,"). ",z$title," *",z$journal,
    "*. [Original research article](",z$journal_article_url,").")
},character(1))
ref_note <- c("## Scientific references: original research only","",
  "The three core references below report original experiments. Lecture PDFs are writing/search guidance and are not scientific references for the abstract. The abstract's own results come from this project's measurements and analyses.","",
  paste(citations,collapse="\n\n"))
writeLines(ref_note,file.path(out,"primary_references.md"))
writeLines(c(readLines(file.path(out,"abstract_draft.md"),warn=FALSE),"",ref_note),
  file.path(out,"abstract_draft.md"))
p <- r$provenance
write_tsv(tibble(source="PubMed MEDLINE metadata",url=p$request_url,
  retrieved_utc=p$retrieved_utc,status=p$http_status,cache_hit=p$cache_hit),file.path(out,"bibliography_API_provenance.tsv"))
message("Draft checked: body=",body_count," words; with title=",body_count+title_count,
  ". Five exact PubMed/DOI records saved for reference-manager import.")
