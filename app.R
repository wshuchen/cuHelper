## A shiny app to help in clinical curation of genomic variants.
## Wenshu Chen, 2026

library(shiny)
library(bslib)
library(shinyjs)
library(grantham)
library(DT)

ui = page_navbar(
              shinyjs::useShinyjs(),
              theme = bs_theme(version = 5, bootswatch = "united"),
              title = h2(HTML("<b>CuHelper</b>"), 
                         style = "font-style: italic; color: green; margin-bottom: 50px;"),

    ## Small tools
    nav_panel(h4(HTML("<b>Tools</b>")),
              
              ## SpliceAI output
              h5(HTML("<b>SpliceAI composer</b>"), 
                 a(id = "toggleSplice", h5(HTML("<b>show/hide</b>")), href = "#")),
              wellPanel(
                  id = "splice",
                  h6("Write SpliceAI search result into one line."),
                  card(
                      layout_columns(
                          col_widths = c(5, 7),
                          card(card_header("SpliceAI scores"),
                               tags$style(type = "text/css", "textarea {width:100%}"),
                               tags$textarea(id = "spliceAIscore",
                                             value = "",
                                             placeholder = "Copy and paste SpliceAI result here.",
                                             rows = 10)
                          ),
                          card(card_header("Result in one line"),
                               uiOutput("spliceAI")
                          )
                      )
                  )),
              
              ## Allele variant nomenclature composer
              h5(HTML("<b>Complex allele composer</b>"), 
                 a(id = "toggleAllele", h5(HTML("<b>show/hide</b>")), href = "#")),
              shinyjs::hidden(wellPanel(
                  id = "allele",
                  h6(markdown("More information on HGVS (https://hgvs-nomenclature.org/stable/)")),
                  layout_columns(
                      col_widths = c(3, 2, 2, 5),
                      card(
                          card_header("variant"),
                          textInput("variant1", "1", 
                                    value = "",
                                    placeholder = "c.1836+1A>G (p.?)"),
                          textInput("variant2", "2", 
                                    value = "",
                                    placeholder = "c.2351C>G (p.Pro784Arg)")
                      ),
                      card(
                          card_header("phase"),
                          radioButtons("phase", "",
                                       choices = list("unknown" = "unknown",
                                                      "in cis" = "cis",
                                                      "in trans" = "trans"),
                                       selected = "unknown")
                      ),
                      card(
                          card_header("protein change"),
                          radioButtons("protein", "",
                                       choices = list("predicted" = "predicted",
                                                      "confirmed" = "confirmed"),
                                       selected = "predicted")
                      ),
                      card(
                          card_header("nomenclature"),
                          uiOutput("nomen")
                      )
                  )
              )),
              
              ## Check if two papers have overlapping authors
              h5(HTML("<b>Check overlapping authors</b>"), 
                 a(id = "toggleAuthor", h5(HTML("<b>show/hide</b>")), href = "#")),
              shinyjs::hidden(wellPanel(
                  id = "author",
                  h6("Provide two PubMed IDs OR two lists of authors. The lists should not include any special symbol."),
                  layout_columns(
                      col_widths = c(3, 6, 3),
                      card(
                          card_header("PubMed ID"),
                          numericInput("id1", "PubMed ID 1",
                                       value = ""),
                          numericInput("id2", "PubMed ID 2",
                                       value = "")
                      ),
                      card(card_header("PubMed authors"),
                           card(
                               tags$style(type = "text/css", "textarea {width:100%}"),
                               tags$textarea(id = "authors1",
                                             value = "",
                                             placeholder = "Copy and paste paper 1 authors here.",
                                             rows = 5)
                           ),
                           card(
                               tags$style(type = "text/css", "textarea {width:100%}"),
                               tags$textarea(id = "authors2",
                                             value = "",
                                             placeholder = "Copy and paste paper 2 authors here.",
                                             rows = 5)
                           ),
                      ),
                      card(card_header("Overlapping authors"),
                           uiOutput("sameauthors")
                      )
                  )
                  
              )),
              
              ## Calculator
              h5(HTML("<b>Simple calculator</b>"), 
                 a(id = "toggleCalculator", h5(HTML("<b>show/hide</b>")), href = "#")),
              shinyjs::hidden(wellPanel(
                  id = "calculator",
                  layout_columns(
                      col_widths = c(5, 2, 5),
                      card(
                          numericInput("number1", "Number1",
                                       value = ""),
                          numericInput("number2", "Number2",
                                       value = "")
                      ),
                      card(
                          radioButtons("operater", "I want to",
                                       choices = list("+" = "+",
                                                      "-" = "-",
                                                      "x" = "*",
                                                      "/" = "/"),
                                       selected = "-")
                      ),
                      card(card_header("Result"),
                           textOutput("result")
                      )
                  )
                  
              ))      
    ),

    ## HGMD ref
    nav_panel(h4(HTML("<b>HGMD</b>")),
              h5(HTML("<b>Convert HGMD reference list to a simpler table.</b>")),
              card(
                  tags$style(type = "text/css", "textarea {width:100%}"),
                  tags$textarea(id = "original",
                                value = "",
                                placeholder = "Copy and paste references here. References should have PubMed ID.",
                                rows = 10)
              ),
              card(
                  layout_columns(
                      col_widths = c(3, 9),
                      card(
                          card(uiOutput("pmid1")),
                          card(uiOutput("pmid2")),
                      ),
                      card(
                          DTOutput("ref_table")
                      )
                  )
              )
    ),
    
    ## Grantham                
    nav_panel(h4(HTML("<b>Grantham</b>")),
              h5(HTML("<b>Amino acid Grantham distance lookup</b>")),
              wellPanel(
                  checkboxGroupInput("aa_x", "Amino acid - X",
                                     choices = sort(c(amino_acids())),
                                     inline = TRUE,
                                     selected = "No"), 
                  checkboxGroupInput("aa_y", "Amino acid - Y",
                                     choices = sort(c(amino_acids())),
                                     inline = TRUE,
                                     selected = "No")
              ),
              wellPanel(
                  textOutput("grantham_distance")
              ),
              wellPanel(
                  DTOutput("aa_property")
              )
    ),
    
    ## Titin
    nav_panel(h4(HTML("<b>Titin</b>")),
              h5(HTML("<b>TTN exon lookup</b>")),
              h6(markdown("Provide _single_ SNV position in hg19 OR hg38, 
                          OR exon range for del-dup.")),
              wellPanel(
                  layout_columns(
                      col_widths = c(6, 6),
                      card(
                          card_header(markdown("**SNV**")),
                          radioButtons("genome", "Genome",
                                       choices = c("H19", "H38"),
                                       selected = "H19"),
                          
                          numericInput("cnv_pos", "Position",
                                       value = ""),
                      ),
                      card(
                          card_header(markdown("**Exon del-dup**")),
                          numericInput("from_exon", "From exon",
                                       value = ""),
                          numericInput("to_exon", "To exon",
                                       value = "")
                      )
                  )
              ),
              wellPanel(
                  h6(markdown("Table adapted from https://www.cardiodb.org/titin/titin_transcripts.php. 
                              Meta to Nvx3 are differant transcripts.")),
                  DTOutput("TTN_exon")
              )
    )
)

server <- function(input, output, session) {
    library(stringr)

    # Grantham distance calculation
    grantham_dist = reactive({
                    D = data.frame()
                    x = c(input$aa_x)
                    y = c(input$aa_y)
                    if (length(x) > 0 & length(y) > 0) {
                        aa_pairs = amino_acid_pairs(x, y)
                        D = grantham_distance(aa_pairs$x, aa_pairs$y)
                    }
                    D
    })
    
    ## Amino acid property display
    aa_df = read.table("data/amino_acid_chemical_property.tsv", 
                       header = TRUE, sep = "\t")
    
    aa_table = reactive({
                AA = data.frame()
                x = c(input$aa_x)
                y = c(input$aa_y)
                if (length(x) > 0) AA = aa_df[aa_df$abbr %in% x, ]
                if (length(y) >0) AA = aa_df[aa_df$abbr %in% y, ]
                if (length(x) > 0 & length(y) >0) AA = aa_df[aa_df$abbr %in% c(x, y), ]
                AA
    })
    
    ## HGMD reference conversion to a table. Four lines for each reference.
    library(rentrez)
    
    ref_table = reactive({
                refs = strsplit(input$original, "\n")[[1]]
                refs = refs[refs != ""]
                if (length(refs) < 4) {
                    ref_df = data.frame(matrix(c(0), ncol = 7))
                    colnames(ref_df) = c("author", "pmid", "year", "disease",
                                         "comment", "PubMed", "journal")
                } else {
                    records = grep("PubMed", refs, value = TRUE)
                    diseases = grep("Disease", refs, value = TRUE)
                    comments = grep("Comments", refs, value = TRUE)
                    pmid = unname(sapply(records, function(x) {
                                        str_trim(str_extract(x, "\\s[0-9]+$"))}))
                    diseases = diseases[!is.na(pmid)]
                    comments = comments[!is.na(pmid)]
                    pmid = pmid[!is.na(pmid)]
                    ref_disease = unname(sapply(diseases, function(x) {
                        str_trim(str_split_i(x, "Disease:", 2))}))
                    ref_comment = unname(sapply(comments, function(x) {
                               str_trim(str_split_i(x, "Comments:", 2))}))
                    
                    # Retrieve info from pubmed summary: author, year, doi
                    esum = entrez_summary(db = "pubmed", id = pmid)
                    if (length(esum$uid) == 1) {
                        pinfo = c(
                                  esum$authors[1, 1],
                                  esum$uid,
                                  substr(esum$pubdate, 1, 4),
                                  esum$articleids[esum$articleids$idtype == "doi", "value"]
                                  )
                        ref_df = data.frame(t(as.data.frame(pinfo)), row.names = esum$uid)
                    } else {
                        pinfo = sapply(esum, function(x) c(
                                     x$authors[1, 1],
                                     x$uid,
                                     substr(x$pubdate, 1, 4),
                                     x$articleids[x$articleids$idtype == "doi", "value"]
                                     )
                                )
                        ref_df = data.frame(t(as.data.frame(pinfo)))
                    }
                    colnames(ref_df) = c("author", "pmid", "year", "doi")
                    
                    pm_url = paste0("https://pubmed.ncbi.nlm.nih.gov/", ref_df$pmid, "/")
                    pm_link = paste0('<a href=', pm_url, ' target="_blank"', '>open</a>')
                    
                    doi_url = paste0("https://doi.org/", ref_df$doi)
                    doi_link = paste0('<a href=', doi_url, ' target="_blank"', '>open</a>')
                    
                    ref_df$disease = ref_disease
                    ref_df$comment = ref_comment
                    ref_df$PubMed = pm_link
                    ref_df$journal = doi_link
                    
                    ref_df$doi = NULL
                    ref_df = ref_df[!duplicated(ref_df), ]
                    ref_df = ref_df[order(ref_df$year), ]
                    rownames(ref_df) = 1:nrow(ref_df)
                }
                ref_df
    })
    
    ## Complex allele variant composer
    cp_nomen = function(variant) {
            # "c.2351C>G (p.Pro784Arg)"
            v = str_trim(str_split_i(variant, "\\(", 1))
            v = str_split_i(v, "c\\.", 2)
            p = str_trim(str_split_i(variant, "p\\.", 2))
            p = str_split_i(p, "\\)", 1)
            return(c(v, p))
    }
    
    write_nomen = function(variant1 = input$variant1, 
                           variant2 = input$variant2, 
                           phase = "unknown", 
                           protein = "predicted") {
        
            v1 = cp_nomen(variant1)
            c1 = v1[1]
            p1 = v1[2]
            v2 = cp_nomen(variant2)
            c2 = v2[1]
            p2 = v2[2]
    
            if (protein == "confirmed") {
                if (phase == "cis") {
                    nomen = paste0("c.[", c1, ";", c2, "] ", "p.[", p1, ";", p2, "]")
                }
                if (phase == "trans") {
                    nomen = paste0("c.[", c1, "];[", c2, "] ", "p.[", p1, "];[", p2, "]")
                }
                if (phase == "unknown")  {
                    nomen = paste0("c.", c1, "(;)", c2, " ", "p.", p1, "(;)", p2)
                }}
            else {
                if (phase == "cis") {
                    nomen = paste0("c.[", c1, ";", c2, "] ", "p.[(", p1, ";", p2, ")]")
                }
                if (phase == "trans") {
                    nomen = paste0("c.[", c1, "];[", c2, "] ", "p.[(", p1, ")];[(", p2, ")]")
                }
                if (phase == "unknown") {
                    nomen = paste0("c.", c1, "(;)", c2, " ", "p.(", p1, ")(;)(", p2, ")")
                }
            }
            return(nomen)
    }
    
    ## Chech overlapping authors
    check_authors_id = function(id1, id2) {
                esum = entrez_summary("pubmed", c(id1, id2))
                same_authors = intersect(esum[[1]]$authors$name, 
                                         esum[[2]]$authors$name)
                return(same_authors)
    }
    
    make_author_list = function(pubmed_authors) {
                authors = gsub("\\d", "", pubmed_authors)
                authors = gsub("\\.", "", pubmed_authors)
                authors = strsplit(authors, ",")
                authors = sapply(authors, function(x) str_trim(x))
                return(authors)
    }

    same_authors = reactive({
                same = ""
                if (!is.na(input$id1) && !is.na(input$id2)) {
                    same = check_authors_id(input$id1, input$id2)
                }
                else {
                    list1 = make_author_list(input$authors1)
                    list2 = make_author_list(input$authors2)
                    if (length(list1) > 1 && length(list2) > 1) {
                        same = intersect(list1,list2)
                    }
                }
                same
    })
    
    ## TTN exon lookup
    ttn_df = read.table("data/ttn.txt", header = TRUE, sep = "\t")
    exon_start = numeric()
    exon_end = numeric()
    exon_df = data.frame()

    ttn_table = reactive({
            if (input$genome == "H19") {
                exon_start = ttn_df$hg19_start
                exon_end = ttn_df$hg19_end
            } else {
                exon_start = ttn_df$hg38_start
                exon_end = ttn_df$hg38_end
            }
            ### CNV
            cnv_pos = input$cnv_pos
            from_exon = input$from_exon
            to_exon = input$to_exon
    
            if (!is.na(cnv_pos)) {
                exon_n = cnv_pos <= exon_start & cnv_pos >= exon_end
                exon_df = ttn_df[exon_n, ]
            }
            if (is.na(cnv_pos) & (!is.na(from_exon) & !is.na(to_exon))) {
                exon_df = ttn_df[from_exon:to_exon, ]  
            }
            exon_df
    })
    
    ## SpliceAI output converter
    splice_scores = reactive({
            if (!nzchar(input$spliceAIscore)) {
                "Please provide SpliceAI search result"
            } else {
                ai_scores = strsplit(input$spliceAIscore, "\n")[[1]]
                scores = lapply(ai_scores, function(x) strsplit(x, "\t")[[1]])
                scores = sapply(scores, function(x) {
                                c(tolower(x[1]), gsub(" ", "", x[3]), x[2])
                })
                scores = data.frame(t(as.data.frame(scores)))
                colnames(scores) = c("loss_gain", "distance", "score")
                scores = scores[!is.na(scores$distance) & 
                                    scores$score > 0 & scores$score != "0.00", ]
                scores = scores[order(scores$score, decreasing = TRUE), ]
                scores = apply(scores, 1, function(x) 
                                c(paste0(x[1], " at ", x[2], ": ", x[3])))
                scores
            }
    })
    
    ## Output
    ## HGMD reference conversion
    observe({
            RT = ref_table()
            if (RT$pmid[1] == 0) {
                output$pmid1 = renderText("PubMed ID output")
            } else {
                output$pmid1 = renderUI({
                    HTML(sapply(RT$pmid, function(x) {
                        paste0("PubMed: ", x, "<br>", sep = "")
                    }))
                })
                output$pmid2 = renderUI({
                    pmids = paste0(RT$pmid, collapse = ", ")
                    HTML(paste0("PubMed: ", pmids))
                    })
            }
            output$ref_table = renderDT(datatable(RT, escape = FALSE))
    })
    
    ## Grantham distance
    observe({
            if (nrow(grantham_dist()) == 0) {
                output$grantham_distance = renderText(
                                        "Please select amino acid pairs."
                                        )
            } else {
                D = grantham_dist()
                output$grantham_distance = renderText(
                    paste(D$x, ">", D$y, "=", D$d, collapse = " | "))
            }
    })
    
    ## Grantham distance matrix
    observe({
            AT = aa_table()
            if (nrow(AT) == 0) {
                output$aa_property = renderDT(aa_df)
            } else {
                output$aa_property = renderDT(AT)
            }
    })
    
    ## TTN
    observe({
            ttn_table = ttn_table()
            if (nrow(ttn_table) >= 1) {
                output$TTN_exon = renderDT(ttn_table)
            } else {
                output$TTN_exon = renderDT(ttn_df)            
            }
    })
    
    ## spliceAI output
    observe({
            splice_scores = splice_scores()
            not_yet = grepl("search", strsplit(splice_scores, " "))
            if (sum(not_yet) == 1) {
                output$spliceAI = renderUI({HTML(paste0(splice_scores, "."))})
            } else {
                splice_scores = paste(splice_scores, collapse = ", ")
                output$spliceAI = renderUI({HTML(paste0("SpliceAI ", splice_scores, "."))})
            }
    })
    
    ## Allele variants
    nomen = reactive({nomen = write_nomen(input$variant1, input$variant2, 
                                          input$phase, input$protein)
            })
    observe({
            nomen = nomen()
            output$nomen = renderUI({HTML(paste(nomen))})
    })
    
    ## overlapping authors
    observe({
            same_authors = same_authors()
            n_same = length(same_authors)
            same_authors = HTML(paste(same_authors, "<br>"))
            if (n_same > 1) {
                output$sameauthors = renderUI({
                    HTML(paste0(n_same, "<br>", same_authors))
                })
            } else {
                output$sameauthors = renderUI({same_authors})            
            }
    })
    
    ## Calculator
    calculator_result = reactive({
            if (is.na(input$number1) && is.na(input$number2)) {""}
            else if (input$operater == "+") {input$number1 + input$number2}
            else if (input$operater == "-") {input$number1 - input$number2}
            else if (input$operater == "*") {input$number1 * input$number2}
            else if (input$operater == "/") {input$number1 / input$number2}
    })
    output$result = renderText(calculator_result())
    
    ## show/hide
    observe({shinyjs::onclick("toggleSplice", 
            shinyjs::toggle(id = "splice", anim = TRUE))
            })
    observe({shinyjs::onclick("toggleAllele", 
            shinyjs::toggle(id = "allele", anim = TRUE))
    })
    observe({shinyjs::onclick("toggleAuthor", 
            shinyjs::toggle(id = "author", anim = TRUE))
    })
    observe({shinyjs::onclick("toggleCalculator", 
            shinyjs::toggle(id = "calculator", anim = TRUE))
    })
}

## Run
shinyApp(ui = ui, server = server)
