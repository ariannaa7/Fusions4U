# Prep - load packages ####
library(shiny) # v1.10.0
library(shinyWidgets) # v0.8.7
library(tidyverse) # v2.0.0
library(DT) # v0.33, for displaying df in multiple pages
library(reactable) # v0.4.4
library(readxl) # v1.4.3

# Error messages
options(shiny.error = traceback)
options(error = function() { traceback(2); quit(save = "no", status = 1, runLast = FALSE) })

# Load & updated necessary datasets ####

# Load complete dataset of validated & annotated fusions
val_ann_fusions_dat <-
  read.csv("fusions_valAnnotated.tsv",
    sep = "\t",
    header = TRUE
  )

# Load the supplementary table to the "cell line info" section
supp_cellLine <-
  read.csv("supplementary_cellLine_info.tsv",
    sep = "\t",
    header = TRUE
  ) %>%
  rename(
    "validated/predicted (arriba)" = "validated.predicted..arriba.",
    "validated/predicted (sf)" = "validated.predicted..sf.",
    "validated/predicted (arriba;sf)" = "validated.predicted..arriba.sf.",
    "discBreakSupport/validated (arriba)" = "discBreakSupport.validated..arriba.",
    "discBreakSupport/validated (sf)" = "discBreakSupport.validated..sf.",
    "discBreakSupport/validated (arriba;sf)" = "discBreakSupport.validated..arriba.sf.",
    "discSupportOnly/validated (arriba)" = "discSupportOnly.validated..arriba.",
    "discSupportOnly/validated (sf)" = "discSupportOnly.validated..sf.",
    "discSupportOnly/validated (arriba;sf)" = "discSupportOnly.validated..arriba.sf."
  )

# Load supplementary excels for "about" section
col_guide <-
  read_excel("column_guide.xlsx", sheet = "colDescriptions", col_names = TRUE)

col_moreInfo <-
  read_excel("column_guide.xlsx", sheet = "moreInfo", col_names = TRUE)

# Load the cell line and gene links files
cell_line_links <-
  read.csv("cell_line_links.tsv",
           sep = "\t",
           header = TRUE
  )

gene_card_links <-
  read.csv("gene_card_links.tsv",
           sep = "\t",
           header = TRUE
  )

# Update the primary df to be shiny ready
shiny_valAnnFusions_df <- val_ann_fusions_dat %>%
  # Join with the gene card links and cell line links df
  left_join(gene_card_links) %>%
  left_join(cell_line_links) %>%
  # Relocate the links columns to be next to the non-link version
  relocate(fusion_name_link, .after = fusion_name) %>%
  relocate(cell_line_link, .after = sampleID) %>%
  arrange(fiveprime_gene_name, threeprime_gene_name, sampleID, disease) # sort the df by fiveprime_gene_name, threeprime_gene_name, sampleID, and disease! Essentially sorting by fusion_name, but avoiding that because of the links in the col now
  
# Identify twin IDs ####
indiv_twin_fusionIDs <- shiny_valAnnFusions_df %>%
  filter(!(is.na(twinID))) %>% # pull the 4788 rows with a twinID ()
  filter(!(str_detect(fusion_id, "_"))) %>% # pull the 3192 rows where twin fusions are represented by two rows
  select(fusion_id)

combo_twin_ids <- shiny_valAnnFusions_df %>%
  filter(!(is.na(twinID))) %>% # pull the 4788 rows with a twinID
  filter((str_detect(fusion_id, "_"))) %>% # pull the 1596 rows where twin fusions are represented by one row
  select(fusion_id)

# Identify genes involved in validated fusions ####
val_gene_names <- shiny_valAnnFusions_df %>%
  filter(program != "both") %>%
  select(fiveprime_gene_name) %>% # 10997, as expected
  rename(gene_name = "fiveprime_gene_name") %>%
  bind_rows(
    shiny_valAnnFusions_df %>% filter(program != "both") %>% select(threeprime_gene_name) %>% dplyr::rename(gene_name = "threeprime_gene_name")
  ) %>% # 21994 as expected (10997 + 10997)
  distinct() # 5976 distinct genes

# Identify recurrent fusions ####
recurrent_fusion_pairs <- shiny_valAnnFusions_df %>%
  filter(program != "both") %>% # 10997, as expected
  select(geneID_fusionName, sampleID) %>%
  distinct() %>% # drops to 5384
  group_by(geneID_fusionName) %>% 
  filter(n() > 1) %>% # 496! regardless of program, identify fusion pairs which appear in more than one cell line!
  select(geneID_fusionName) %>%
  distinct() %>% # 185 fusion pairs which appear accross multiple different cell lines
  pull(geneID_fusionName) # now a vector of the geneID_fusionNames we want

# Create a map type df with the tissue type as it appears in the table and a column for how we want it displayed in the drop down
tissue_aliases <- shiny_valAnnFusions_df %>%
  select(tissueType) %>%
  distinct() %>%
  arrange(tissueType) %>%
  mutate(tissue_alias = case_when(
    (tissueType == "AUTONOMIC_GANGLIA") ~ "Autonomic ganglia",
    (tissueType == "BONE") ~ "Bone",
    (tissueType == "BREAST") ~ "Breast",
    (tissueType == "CENTRAL_NERVOUS_SYSTEM") ~ "Central nervous system",
    (tissueType == "ENDOMETRIUM") ~ "Endometrium",
    (tissueType == "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE") ~ "Haematopoietic & lymphoid tissue",
    (tissueType == "KIDNEY") ~ "Kidney",
    (tissueType == "LARGE_INTESTINE") ~ "Large intestine",
    (tissueType == "LIVER") ~ "Liver",
    (tissueType == "LUNG") ~ "Lung",
    (tissueType == "OESOPHAGUS") ~ "Oesophagus",
    (tissueType == "OVARY") ~ "Ovary",
    (tissueType == "PANCREAS") ~ "Pancreas",
    (tissueType == "PLEURA") ~ "Pleura",
    (tissueType == "PROSTATE") ~ "Prostate",
    (tissueType == "SKIN") ~ "Skin",
    (tissueType == "SMALL_INTESTINE") ~ "Small intestine",
    (tissueType == "SOFT_TISSUE") ~ "Soft tissue",
    (tissueType == "STOMACH") ~ "Stomach",
    (tissueType == "THYROID") ~ "Thyroid",
    (tissueType == "UPPER_AERODIGESTIVE_TRACT") ~ "Upper aerodigestive tract",
    (tissueType == "URINARY_TRACT") ~ "Urinary tract"
  )) %>%
  rename(tissue_original = "tissueType")

# Create a map type df with the disease as it appears in the table and a column for how we want it displayed in the drop down
disease_aliases <- shiny_valAnnFusions_df %>%
  select(disease) %>%
  distinct() %>%
  arrange(disease) %>%
  mutate(disease_alias = case_when(
    (disease == "carcinoma") ~ "carcinoma",
    (disease == "carcinoma (adenocarcinoma)") ~ "adenocarcinoma",
    (disease == "carcinoma (anaplastic_carcinoma)") ~ "anaplastic carcinoma",
    (disease == "carcinoma (bronchioloalveolar_adenocarcinoma)") ~ "bronchioloalveolar adenocarcinoma",
    (disease == "carcinoma (clear_cell_carcinoma)") ~ "clear cell carcinoma",
    (disease == "carcinoma (clear_cell_renal_cell_carcinoma)") ~ "clear cell renal cell carcinoma",
    (disease == "carcinoma (diffuse_adenocarcinoma)") ~ "diffuse adenocarcinoma",
    (disease == "carcinoma (ductal_carcinoma)") ~ "ductal carcinoma",
    (disease == "carcinoma (ductal_carcinoma\\, medullary)") ~ "medullary ductal carcinoma",
    (disease == "carcinoma (endometrioid_carcinoma)") ~ "endometrioid carcinoma",
    (disease == "carcinoma (hepatocellular_carcinoma)") ~ "hepatocellular carcinoma",
    (disease == "carcinoma (intestinal_adenocarcinoma)") ~ "intestinal adenocarcinoma",
    (disease == "carcinoma (large_cell_carcinoma)") ~ "large cell carcinoma",
    (disease == "carcinoma (metaplastic_carcinoma)") ~ "metaplastic carcinoma",
    (disease == "carcinoma (mixed_adenosquamous_carcinoma)") ~ "mixed adenosquamous carcinoma",
    (disease == "carcinoma (mixed_carcinoma)") ~ "mixed carcinoma",
    (disease == "carcinoma (mucinous_carcinoma)") ~ "mucinous carcinoma",
    (disease == "carcinoma (non_small_cell_carcinoma)") ~ "non-small cell carcinoma",
    (disease == "carcinoma (renal_cell_carcinoma)") ~ "renal cell carcinoma",
    (disease == "carcinoma (serous_carcinoma)") ~ "serous carcinoma",
    (disease == "carcinoma (signet_ring_adenocarcinoma)") ~ "signet ring adenocarcinoma",
    (disease == "carcinoma (small_cell_carcinoma)") ~ "small cell carcinoma",
    (disease == "carcinoma (squamous_cell_carcinoma)") ~ "squamous cell carcinoma",
    (disease == "carcinoma (transitional_cell_carcinoma)") ~ "transitional cell carcinoma",
    (disease == "carcinoma (tubular_adenocarcinoma)") ~ "tubular adenocarcinoma",
    (disease == "carcinoma (undifferentiated_adenocarcinoma)") ~ "undifferentiated adenocarcinoma",
    (disease == "carcinoma (undifferentiated_carcinoma)") ~ "undifferentiated carcinoma",
    (disease == "glioma") ~ "glioma",
    (disease == "glioma (astrocytoma)") ~ "astrocytoma",
    (disease == "glioma (astrocytoma_Grade_III\\, anaplastic)") ~ "anaplastic astrocytoma grade III",
    (disease == "glioma (astrocytoma_Grade_IV\\, glioblastoma_multiforme)") ~ "astrocytoma grade IV",
    (disease == "glioma (gliosarcoma)") ~ "gliosarcoma",
    (disease == "haematopoietic_neoplasm (acute_myeloid_leukaemia\\, M5)") ~ "acute myeloid leukaemia M5",
    (disease == "haematopoietic_neoplasm (acute_myeloid_leukaemia\\, M6)") ~ "acute myeloid leukaemia M6",
    (disease == "haematopoietic_neoplasm (blast_phase_chronic_myeloid_leukaemia)") ~ "blast phase chronic myeloid leukaemia",
    (disease == "leiomyosarcoma") ~ "leiomyosarcoma",
    (disease == "lymphoid_neoplasm (anaplastic_large_cell_lymphoma)") ~ "anaplastic large cell lymphoma",
    (disease == "lymphoid_neoplasm (diffuse_large_B_cell_lymphoma)") ~ "diffuse large B-cell lymphoma",
    (disease == "lymphoid_neoplasm (mantle_cell_lymphoma)") ~ "mantle cell lymphoma",
    (disease == "lymphoid_neoplasm (plasma_cell_myeloma)") ~ "plasma cell myeloma",
    (disease == "malignant_melanoma") ~ "malignant melanoma",
    (disease == "mesothelioma") ~ "mesothelioma",
    (disease == "neuroblastoma") ~ "neuroblastoma",
    (disease == "osteosarcoma") ~ "osteosarcoma",
    (disease == "primitive_neuroectodermal_tumour-medulloblastoma") ~ "medulloblastoma",
    (disease == "rhabdoid_tumour") ~ "rhabdoid tumour",
    (disease == "rhabdomyosarcoma") ~ "rhabdomyosarcoma",
    (disease == "sarcoma") ~ "sarcoma"
  )) %>%
  rename(disease_original = "disease")



# Create vectors with column names which don't need to displayed in shiny app & those that don't need to be downloaded ####
# These cols don't need to be displayed because they are irrelevant (since their info is already included somewhere else) or they are of niche interest and take up too much space in a cell to be included here but can be found when downlaoding the df from shiny
doNot_displayCols <-
  c(
    "fusion_id",
    "fusion_name",
    "sampleID",
    "rna_run_accession",
    "wgs_run_accession",
    "multi_alt_fusion_allPred",
    "altSplicing_id_allPred",
    "fiveprime_gene_name",
    "fiveprime_geneID",
    "fiveprime_gene_strand",
    "fiveprime_transcribed_strand",
    "hg38_fiveprime_chr",
    "hg38_fiveprime_gene_gtfStart",
    "hg38_fiveprime_fusion_junction",
    "hg38_fiveprime_gene_gtfEnd",
    "hg38_fiveprime_searchStart",
    "hg38_fiveprime_searchEnd",
    "hg19_fiveprime_chr",
    "hg19_fiveprime_geneStart",
    "hg19_fiveprime_fusion_junction",
    "hg19_fiveprime_geneEnd",
    "hg19_fiveprime_searchStart",
    "hg19_fiveprime_searchEnd",
    "threeprime_gene_name",
    "threeprime_geneID",
    "threeprime_gene_strand",
    "threeprime_transcribed_strand",
    "hg38_threeprime_chr",
    "hg38_threeprime_gene_gtfStart",
    "hg38_threeprime_fusion_junction",
    "hg38_threeprime_gene_gtfEnd",
    "hg38_threeprime_searchStart",
    "hg38_threeprime_searchEnd",
    "hg19_threeprime_chr",
    "hg19_threeprime_geneStart",
    "hg19_threeprime_fusion_junction",
    "hg19_threeprime_geneEnd",
    "hg19_threeprime_searchStart",
    "hg19_threeprime_searchEnd",
    "hg19_fiveprime_gene_coordinates",
    "hg19_fiveprime_junction_coordinate",
    "hg19_threeprime_gene_coordinates",
    "hg19_threeprime_junction_coordinate",
    "hg19_fiveprime_breakpoint_sanityPassingPair",
    "hg19_threeprime_breakpoint_sanityPassingPair",
    "hg19_fiveprime_breakpoint_sanityFailingPair",
    "hg19_threeprime_breakpoint_sanityFailingPair",
    "hg38_fiveprime_breakpoint_sanityFailingPair",
    "hg38_threeprime_breakpoint_sanityFailingPair",
    "miRNA_fiveprime_gffStart",
    "miRNA_fiveprime_gffEnd",
    "miRNA_threeprime_gffStart",
    "miRNA_threeprime_gffEnd",
    "miRNA_downstream_threeprime_gffStart",
    "miRNA_downstream_threeprime_gffEnd",
    "retained_protein_domains_PfamArriba_arriba",
    "retained_protein_domains_PfamUCSC_arriba", # these pfam/annots columns can make rows width very large col heights. User can still see these upon download
    "fusion_transcript_arriba",
    "peptide_sequence_arriba",
    "read_identifiers_arriba",
    "fiveprime_direction_arriba",
    "threeprime_direction_arriba",
    "filters_arriba",
    "fiveprime_transcript_id_arriba",
    "threeprime_transcript_id_arriba",
    "tags_arriba",
    "spliceType_sf",
    "largeAnchorSupport_sf",
    "annots_sf",
    "fiveprime_BreakDinuc_sf",
    "fiveprime_BreakEntropy_sf",
    "threeprime_BreakDinuc_sf",
    "threeprime_BreakEntropy_sf",
    "fiveprime_cds_ID_sf",
    "fiveprime_cds_range_sf",
    "threeprime_cds_ID_sf",
    "threeprime_cds_range_sf",
    "est_J_sf",
    "est_S_sf",
    "fiveprime_pfam_sf",
    "threeprime_pfam_sf",
    "fusion_model_sf",
    "fusion_cds_sf",
    "fusion_transl_sf",
    "tf_hg19_threeprime_fusion_junction",
    "tf_hg19_fiveprime_fusion_junction",
    "tcga_hg38_threeprime_fusion_junction",
    "tcga_hg38_fiveprime_fusion_junction",
    "tf_threeprime_exactOrNear_juncMatch",
    "tf_fiveprime_exactOrNear_juncMatch",
    "tcga_fiveprime_exactOrNear_juncMatch",
    "tcga_threeprime_exactOrNear_juncMatch",
    "mitelman_tissue_match",
    "tf_tissue_match",
    "tcga_tissue_match",
    "tcga_tissueType",
    "tf_tissueType"
  )

doNotDisplayCols_butNeedInDetailView <-
  c(
    "geneID_fusionName",
    "disease_stage",
    "twinID",
    "multi_alt_fusion_valOnly",
    "altSplicing_id_onlyVal",
    "fiveprime_gene_type",
    "hg38_fiveprime_gene_gtf_coordinates",
    "threeprime_gene_type",
    "hg38_threeprime_gene_gtf_coordinates",
    #"hg38_fiveprime_breakpoint_sanityPassingPair", # keeping these now, adding elipses in table view the case of multiple breakpoints
    #"hg38_threeprime_breakpoint_sanityPassingPair",
    "fiveprime_transcriptionType",
    "threeprime_transcriptionType",
    "fusionType_arriba",
    "discordant_mates_arriba",
    "fiveprime_split_reads_arriba",
    "threeprime_split_reads_arriba",
    "fiveprime_coverage_arriba",
    "threeprime_coverage_arriba",
    "junctionReadCount_sf",
    "spanningFragCount_sf",
    "ffpm_sf",
    "fiveprime_miRNA_host",
    "miRNA_fiveprime_coords",
    "threeprime_miRNA_host",
    "miRNA_threeprime_coords",
    "threeprime_candidate_miRNA_host",
    "miRNA_downstream_threeprime_coords",
    "fiveprime_kinase",
    "fiveprime_inFrame_kinase",
    "fiveprime_kinase_family",
    "fiveprime_kinase_subFamily",
    "fiveprime_kinase_group",
    "fiveprime_kinase_uniprotID",
    "threeprime_kinase",
    "threeprime_inFrame_kinase",
    "threeprime_kinase_family",
    "threeprime_kinase_subFamily",
    "threeprime_kinase_uniprotID",
    "threeprime_kinase_group",
    "mitelmanFusion_match",
    "tf_match",
    "tf_position_consistency",
    "tf_frame_prediction",
    "tcga_match",
    "cgc_fiveprime_gene",
    "cgc_fiveprime_tier",
    "cgc_fiveprime_hallmark",
    "cgc_fiveprime_tissue_type",
    "cgc_fiveprime_tumour_types_somatic",
    "cgc_fiveprime_tumour_types_germline",
    "cgc_threeprime_gene",
    "cgc_threeprime_tier",
    "cgc_threeprime_hallmark",
    "cgc_threeprime_tissue_type",
    "cgc_threeprime_tumour_types_somatic",
    "cgc_threeprime_tumour_types_germline"
  )
  

# These cols don't need to be downloaded (because you can already pull their info from another col)
doNot_downloadCols <-
  c(
    "fusion_name_link",
    "cell_line_link", # don't download the long urls, download the actual fusion name and cell line/sample ID cols
    "hg38_fiveprime_searchStart",
    "hg38_fiveprime_searchEnd",
    "hg19_fiveprime_searchStart",
    "hg19_fiveprime_searchEnd",
    "hg38_threeprime_searchStart",
    "hg38_threeprime_searchEnd",
    "hg19_threeprime_searchStart",
    "hg19_threeprime_searchEnd",
    "miRNA_fiveprime_gffStart",
    "miRNA_fiveprime_gffEnd",
    "miRNA_threeprime_gffStart",
    "miRNA_threeprime_gffEnd",
    "miRNA_downstream_threeprime_gffStart",
    "miRNA_downstream_threeprime_gffEnd"
  )

# User Interface ####
ui <- fluidPage(

  titlePanel(
    "Fusions4U: RNA-Seq predicted and WGS validated gene fusions (hg38)"
  ),
  sidebarLayout(
    sidebarPanel(
      # On default, show all fusions. Allow user to pick 5' gene, 3' gene, or both. When they pick only one of the genes or none of the genes, the default will be "All"
      width = 3, # width of panel
      div(
        style = "display: flex; justify-content: center; gap: 2px;", # put the buttons NEXT to each other instead of one below the other, center the buttons to width of sidebar, and put space between them
        # trial and error with font size, style (above), and white-space/flex/width options. This displays well in smaller laptop screen and larger display
        # User can click this button to clear all side panel filtering and reset to original
        actionButton("reset", "Reset to defaults", class = "btn btn-danger", style = "font-size: 85%; white-space: nowrap; flex: 1;"), # btn-danger makes the button red
        
        # User can click this button to download fusions table (filtered/unfiltered depending on sidebar filtering user has utilized)
        downloadButton("download_fusions_side", "Download fusions (TSV)", class = "btn btn-primary", style = "font-size: 85%; white-space: nowrap; flex: 1;") # btn-primary makes the button blue
      ),

      # Display twins with one row or separate rows for Arriba & STAR-Fusion
      radioButtons(
        "format",
        "Display identical Arriba & STAR-Fusion validated predictions in:",
        choices = c(
          "Separate rows" = "separate",
          "One row" = "combo"
        ),
        selected = "separate"
      ),

      # User can choose to display all fusions or just those identically predicted
      radioButtons(
        "type",
        "Validated fusions to display:",
        choices = c(
          "All" = "all",
          "Only those identically predicted by Arriba & STAR-Fusion" = "twins"
        ),
        selected = "all"
      ),

      # Allow user to select if they would like to specify a gene as a specific partner OR just select a gene of interest
      radioButtons(
        "gene_selection",
        "Gene filtering options:",
        # Drop down for 5' partner and drop down for 3' partner
        choices = c(
          "Select gene partners" = "both",
          "Select a gene of interest" = "one"
        ),
        selected = "both" # just show all genes, could be either a 5' or 3' partner just important that it is involved in the fusion
      ),
    
      # Conditional for gene_selection radio buttons
      conditionalPanel(
        condition = "input.gene_selection == 'both'",
        div(
          style = "margin-left: 20px;", # indent the drop down to be nestled under the associated radio buttons
          # Allow user to select a 5' partner from a dropdown (can also type to narrow results)
          selectizeInput(
            "fivePrime",
            "5' gene:",
            # choices = c("All", unique(
            #   shiny_valAnnFusions_df$fiveprime_gene_name
            # )),
            choices = c("All"),
            multiple = FALSE,
            selected = "All",
            options = list(
              placeholder = "Type to search...",
              allowEmptyOption = TRUE
            )
          ),
          # Allow user to select a 3' partner from a dropdown (can also type to narrow results)
          selectizeInput(
            "threePrime",
            "3' gene:",
            # choices = c("All", unique(
            #   shiny_valAnnFusions_df$threeprime_gene_name
            # )),
            choices = c("All"),
            multiple = FALSE,
            selected = "All",
            options = list(
              placeholder = "Type to search...",
              allowEmptyOption = TRUE
            )
          ) 
        )
      ),
      conditionalPanel(
        condition = "input.gene_selection == 'one'",
        div(
          style = "margin-left: 20px;", # indent the drop down to be nestled under the associated radio buttons
          # allow user to select A partner (as specified, 5' or 3') from a dropdown (can also type to narrow results)
          selectizeInput(
            "gene",
            "Gene:",
            # choices = c("All", val_gene_names$gene_name),
            choices = c("All"),
            multiple = FALSE,
            selected = "All",
            options = list(
              placeholder = "Type to search...",
              allowEmptyOption = TRUE
            )
          ) 
        )
      ),

      # Filter by cell line!
      selectizeInput(
        "cellLine",
        "Select cell line of interest",
        choices = NULL, # empty initially
        multiple = FALSE,
        selected = "All"
      ),
      
      # Filter by tissue type of cell line! Can type to narrow selection and can select more than one tissue
      selectizeInput(
        "tissue",
        "Select tissue(s) of interest",
        choices = NULL, # empty initially
        multiple = TRUE,
        selected = NULL
      ),

      # Filter by disease of cell line! Can type to narrow selection and can select more than one diseases
      selectizeInput(
        "disease",
        "Select disease(s) of interest",
        choices = NULL, # empty initially
        multiple = TRUE,
        selected = NULL
      ),

      # Select radio button to display fusions of certain frame
      # Change to radio! Checkbox allows to select more than one, then looks for one fusion that is both in-frame and out-of-frame (not possible) instead of showing those fusions that are in-frame and those that are out-of-frame.
      radioButtons(
        "readingFrame",
        "Predicted reading frame",
        choices = list(
          'Any (includes "unknown")' = "any",
          "In-frame" = "in_frame",
          "Out-of-frame" = "out_of_frame",
          "Stop codon (Arriba only)" = "stop_codon"
        ),
        selected = "any"
      ),

      # Select checkboxes to display fusions only with certain annotations
      checkboxGroupInput(
        "annotations",
        "Filter by annotation type",
        choices = list(
          "Alternative splicing support" = "alt",
          "Recurring fusion pairs" = "recur",
          "Possible promoter swapping event" = "promoter",
          "Kinase partner" = "kinase",
          "miRNA host gene parter" = "miRNA",
          "Fusion parter present in COSMIC Cancer Gene Census" = "cgc",
          "Mitelman fusion match" = "mitelman",
          "TumorFusions match" = "tf",
          "TCGA - validated FusionCatcher match" = "tcga"
        ),
        selected = NULL
      ),
    ),
    mainPanel(
      uiOutput("instructions_blurb"),
      tabsetPanel(
        # ID for tabs, call it down the line as input$tabs_id
        id = "tabs_ID",
        tabPanel(
          "Table View",
          DTOutput("table"), # table that is searchable, and allows to click for next page
          downloadButton("download_fusions_bottom", "Download fusions (TSV)", class = "btn btn-primary") # button BELOW TABLE that allows user to download the table
        ),
        tabPanel(
          "Detail View", # select a row in "Table View" and see the info more easily here
          fluidRow(
            column(6, uiOutput("topLeft_card")), # general info section
            column(6, uiOutput("topRight_card"))
          ),
          # Additional info section
          fluidRow(
            column(6, uiOutput("bottomLeft_card")), # 5' partner section
            column(6, uiOutput("bottomRight_card")) # 3' partner section
          ) 
        ),
        tabPanel(
          "Cell Line Supplementary", # see supplementary table that was created with rows for each cell line
          uiOutput("supp_blurb"),
          DTOutput("supp_table"), # table that is searchable, and allows to click for next page
          downloadButton("download_supp", "Download supplementary (TSV)", class = "btn btn-primary") # button that allows user to download the supplementary table
        ),
        tabPanel(
          "About", # check here for more information about each column in "Table View"
          uiOutput("col_blurb"),
          uiOutput("more_infoTable"), # table with hyperlinks
          DTOutput("col_table") # table that is searchable, and allows to click for next page
        )
      )
    )
  )
)

# Server ####
server <- function(input, output, session) {
  
  # Clear filters option ####
  observeEvent(input$reset, {# reset to defaults and clear filtering
    updateSelectInput(session, "format", selected = "separate")
    updateRadioButtons(session, "type", selected = "all")
    updateRadioButtons(session, "gene_selection", selected = "both")
    updateSelectInput(session, "fivePrime", selected = "All")
    updateSelectInput(session, "threePrime", selected = "All")
    updateSelectInput(session, "gene", selected = "All")
    updateRadioButtons(session, "cellLine", selected = "All")
    updateCheckboxGroupInput(session, "disease", selected = character(0)) # selected = character(0) resets to NULL
    updateCheckboxGroupInput(session, "tissue", selected = character(0)) # selected = character(0) resets to NULL
    updateCheckboxGroupInput(session, "readingFrame", selected = "any") # selected = character(0) resets to NULL
    updateCheckboxGroupInput(session, "annotations", selected = character(0)) # selected = character(0) resets to NULL
  })
  
  # Instructions blurb ####
  output$instructions_blurb <- renderUI({
    if (input$tabs_ID == "Table View") { # when Table View is selected
      tagList(
        p("Welcome to Fusions4U!"),
        p("Here you will find gene fusions which were predicted using RNA-Seq and further validated with WGS data. Fusions were predicted using Arriba and STAR-Fusion."),
        p("To sort the table, select the arrows beside the column name. To sort by more than one column, hold the SHIFT key while making your selection."),
        p("For a clearer and more informative view, select the row you are interested in and move to the 'Detail View' tab"),
        p("To download the table (filtered or unfiltered), select 'Download fusions (TSV)' in the side panel or at the bottom of your screen. The TSV file will contain many more columns than shown here allowing you to select columns for your interest/use.")
      )
    } else {
      NULL # no blurb for other tabs (i.e. no blurb for Detail View)
    }
  })
  
  # Organize fusion datatable based on user's "format" selection (twins together or sep) ####
  format_df <- reactive({
    df <- shiny_valAnnFusions_df
    
    # Conditionally filter based on the selected option
    if (input$format == "separate") {
      df <- df %>%
        anti_join(combo_twin_ids, by = "fusion_id") # remove the rows with twins rep by one row. This will not affect twins rep by two rows
    }
    
    if (input$format == "combo") {
      df <- df %>%
        anti_join(indiv_twin_fusionIDs, by = "fusion_id") # remove the rows with twins rep on two rows. This will not affect twins rep by one row
    }
    
    df # Return the filtered/reformatted dataset
  })
  
  # Organize fusion datatable based on user's "type" selection (all or only twins) ####
  type_df <- reactive({
    df <- format_df()
    
    # Conditionally filter based on the selected option
    if (input$type == "all") {
      df <- df
    }
    
    if (input$type == "twins") {
      df <- df %>%
        filter(twin == "TRUE") # only keep instances of twins
    }
    
    df # Return the filtered/reformatted dataset
  })
  
  # Filter fusion datatable based on user's gene or fiveprime/threeprime selection ####
  fusionPartner_filtered_df <- reactive({
    df <- type_df()
    
    # If the user chooses not to display all 5' partners, but rather pick one
    if (input$gene_selection == "both" &
        input$fivePrime != "All") {
      df <- df %>%
        filter(fiveprime_gene_name == input$fivePrime) # then filter the df to only show fusions with the partner they picked
    }
    # If the user chooses not to display all 3' partners, but rather pick one
    if (input$gene_selection == "both" &
        input$threePrime != "All") {
      df <- df %>%
        filter(threeprime_gene_name == input$threePrime) # then filter the df to only show fusions with the partner they picked
    }
    # If the user chooses not to display all 3' partners, but rather pick one
    if (input$gene_selection == "one" &
        input$gene != "All") {
      df <- df %>%
        filter(fiveprime_gene_name == input$gene |
                 threeprime_gene_name == input$gene) # then filter the df to only show fusions with the partner they picked
    }
    
    df # just show it all (don't filter) if "All" is selected
  })
  
  # Update threeprime choices based on fiveprime selection, update the type_df() rather than the fusion_filtered one so that we maintain the drop down list even after filtering for genes
  observe({ # updates available options
    if (input$gene_selection != "both") return()
    
    available_threePrime <- if (input$fivePrime == "All") {
      sort(unique(type_df()$threeprime_gene_name))
    } else {
      sort(unique(type_df()$threeprime_gene_name[type_df()$fiveprime_gene_name == input$fivePrime]))
    }
    
    updateSelectizeInput( # updates the drop down
      session,
      "threePrime",
      choices = c("All", available_threePrime),
      selected = if (nzchar(input$threePrime) && input$threePrime %in% available_threePrime) {
        input$threePrime
      } else if (nzchar(input$threePrime)) {
        "All"
      } else {
        "" # empites drop down while typing
      }
    )
      
  })
  
  # Update fiveprime choices based on fiveprime selection, update the type_df() rather than the fusion_filtered one so that we maintain the drop down list even after filtering for genes
  observe({ # updates available options
    if (input$gene_selection != "both") return()
    
    available_fivePrime <- if (input$threePrime == "All") {
      sort(unique(type_df()$fiveprime_gene_name))
    } else {
      sort(unique(type_df()$fiveprime_gene_name[type_df()$threeprime_gene_name == input$threePrime]))
    }
    
    updateSelectizeInput( # updates the dropdown 
      session,
      "fivePrime",
      choices = c("All", available_fivePrime),
      selected = if (nzchar(input$fivePrime) && input$fivePrime %in% available_fivePrime) {
        input$fivePrime
      } else if (nzchar(input$fivePrime)) {
        "All"
      } else {
        "" # empties while typing
      }
    )
  })
  
  # Update the drop down gene choices based on type_df()
  observe({
    if (input$gene_selection != "one") return()
    
    available_genes <- sort(unique(c(type_df()$fiveprime_gene_name, type_df()$threeprime_gene_name)))
    
    updateSelectizeInput(
      session,
      "gene",
      choices = c("All", available_genes),
      selected = if (nzchar(input$gene) && input$gene %in% available_genes) {
        input$gene
      } else if (nzchar(input$gene)) {
        "All"
      } else {
        ""  
      }
    )
  })

  
  # Update selectizeInput for cell line when fusionPartner_filtered_df changes
  observeEvent(fusionPartner_filtered_df(), {
    updateSelectizeInput(
      session,
      "cellLine",
      choices = c("All", sort(unique(fusionPartner_filtered_df()$sampleID))),
      selected = "All"
    )
  })
  
  # Filter based on the selected cell line
  cellLine_filtered_df <- reactive({
    df <- fusionPartner_filtered_df()
    if (input$cellLine == "All") {
      df <- df
    }

    if (input$cellLine != "All") {
      df <- df %>%
        filter(sampleID == input$cellLine) # only keep instances of the selected cell line
    }

    df # Return the filtered/reformatted dataset
  })
  
  # Update selectizeInput for tissue when cellLine_filtered_df changes
  observeEvent(cellLine_filtered_df(), {

    available_tissues <- sort(unique(cellLine_filtered_df()$tissueType)) # tissues in the df
    
    available_tissues_alias <- tissue_aliases[tissue_aliases$tissue_original %in% available_tissues, ] # the associated aliases which we want to display
    
    
    updateSelectizeInput(
      session,
      "tissue",
      choices = sort(unique(available_tissues_alias$tissue_alias)),
      selected = NULL
    )
  })
  
  # Use the tissue mapping df to map from the selected tissue back to how it is written in the table
  selected_tissue_originalName <- reactive({
    req(input$tissue)  # make sure something is selected
    tissue_aliases$tissue_original[
      match(input$tissue, tissue_aliases$tissue_alias)
    ]
  })
  
  # Filter the fusion datatable based on user's "tissue" selection ####
  tissue_filtered_df <- reactive({
    df <- cellLine_filtered_df()
    
    # If nothing selected, return all rows
    if (is.null(input$tissue) || length(input$tissue) == 0) {
      return(df)
    }
    
    # Keep only rows that match the tissue type selected (using the original name)
    df[df$tissueType %in% selected_tissue_originalName(), ]
    
  })
  
  # Update selectizeInput for disease when tissue_filtered_df changes
  observeEvent(tissue_filtered_df(), {
    
    available_diseases <- sort(unique(tissue_filtered_df()$disease)) # tissues in the df
    
    available_diseases_alias <- disease_aliases[disease_aliases$disease_original %in% available_diseases, ] # the associated aliases which we want to display
    
    
    updateSelectizeInput(
      session,
      "disease",
      choices = sort(unique(available_diseases_alias$disease_alias)),
      selected = NULL
    )
  })
  
  # Use the tissue mapping df to map from the selected tissue back to how it is written in the table
  selected_disease_originalName <- reactive({
    req(input$disease)  # make sure something is selected
    disease_aliases$disease_original[
      match(input$disease, disease_aliases$disease_alias)
    ]
  })
  
  # Filter the fusion datatable based on user's "disease" selection ####
  disease_filtered_df <- reactive({
    df <- tissue_filtered_df()
    
    # If nothing selected, return all rows
    if (is.null(input$disease) || length(input$disease) == 0) {
      return(df)
    }
    
    # Keep only rows that match the disease selected (using the original name)
    df[df$disease %in% selected_disease_originalName(), ]
    
  })
  
  
  # Filter the fusion datatable based on user's "readingFrame" selection ####

  frame_filtered_df <- reactive({
    df <- disease_filtered_df()


    if (input$readingFrame == "any") {
      df <- df
    } else if (input$readingFrame == "in_frame") {
      df <- df %>%
        filter(reading_frame_arriba == "in-frame" |
          reading_frame_sf == "in-frame")
    } else if (input$readingFrame == "out_of_frame") {
      df <- df %>%
        filter(reading_frame_arriba == "out-of-frame" |
          reading_frame_sf == "out-of-frame")
    } else if (input$readingFrame == "stop_codon") {
      df <- df %>%
        filter(reading_frame_arriba == "stop-codon")
    }

    df
  })
  
  # Filter the fusion datatable based on user's "annotations" selection ####
  annotation_filtered_df <- reactive({
    df <- frame_filtered_df()

    if ("alt" %in% input$annotations) {
      df <- df %>%
        filter(multi_alt_fusion_valOnly == "TRUE")
    }
    
    if ("recur" %in% input$annotations) {
      df <- df %>%
        filter(geneID_fusionName %in% recurrent_fusion_pairs)
    }
    
    if ("promoter" %in% input$annotations) {
      df <- df %>%
        filter(promoterSwapEventCandidate == "TRUE")
    }

    if ("kinase" %in% input$annotations) {
      df <- df %>%
        filter(fiveprime_kinase == "TRUE" |
                 threeprime_kinase == "TRUE")
    }
    
    if ("miRNA" %in% input$annotations) {
      df <- df %>%
        filter(
          fiveprime_miRNA_host == "TRUE" |
            threeprime_miRNA_host == "TRUE" |
            threeprime_candidate_miRNA_host == "TRUE"
        )
    }
    
    if ("cgc" %in% input$annotations) {
      df <- df %>%
        filter(cgc_fiveprime_gene == "TRUE" |
                 cgc_threeprime_gene == "TRUE")
    }

    if ("mitelman" %in% input$annotations) {
      df <- df %>%
        filter(mitelmanFusion_match == "TRUE")
    }

    if ("tf" %in% input$annotations) {
      df <- df %>%
        filter(tf_match == "TRUE")
    }

    if ("tcga" %in% input$annotations) {
      df <- df %>%
        filter(tcga_match == "TRUE")
    }

    df
  })
  
  # Update the DF based on how it should be displayed
  complete_filtered_df_toDisplay <- reactive({
    df <- annotation_filtered_df() %>%
      select(-all_of(doNot_displayCols)) # omit the columns named in this vector

    df
  })

  # Reformat the filtered fusion datatable to not download certain unnecessary columns ####
  complete_filtered_df_to_download <- reactive({
    df <- annotation_filtered_df() %>%
      select(-all_of(doNot_downloadCols)) # omit the columns named in this vector

    df
  })
  
  # Display the fusions in table format ("Table View") based on filtering/reformatting parameters ####
  output$table <- renderDT({
    
    # Display complete_filtered_df_toDisplay, manually omit geneID_fusionName, and update the "link" style columns to be named like the regular version
    data <- complete_filtered_df_toDisplay() %>% select(-all_of(doNotDisplayCols_butNeedInDetailView)) %>% rename(
      "fusion_name" = "fusion_name_link",
      "sampleID" = "cell_line_link",
      "hg38_fiveprime_breakpoint" = "hg38_fiveprime_breakpoint_sanityPassingPair", # rename for user-friendly purposes. column will maintain original naming in the downloaded version
      "hg38_threeprime_breakpoint" = "hg38_threeprime_breakpoint_sanityPassingPair"
    )
    
    # Create a function to truncate the breakpoint columns after the first comma (if there is one) and add an elpises! Keep function here instead of top of script
    trunc_break_col <- function(breakpoint_col) {
      sapply(breakpoint_col, function(x) {
        if (grepl(",", x)) { # if the cell has a comma, e.g "chr12:53314072,chr12:53314077"
          breakpoints <- strsplit(x, ",")[[1]] # then split the into a vector based on the commas, e.g. c("chr12:53314072","chr12:53314077")
          paste0(
            breakpoints[1], # display the first breakpoint
            ' <span class="expand" style="color:blue; cursor:pointer;">…</span>', # add a clickable ellipsis
            '<span class="full" style="display:none;">, ', paste(breakpoints[-1], collapse = ", "), '</span>' # hide the rest of the breakpoints
          )
        } else {
          x # return original cell if there is no comma in it
        }
      }, USE.NAMES = FALSE)
    }
    
    # Run the function to format the breakpoint cols
    data$hg38_fiveprime_breakpoint <- trunc_break_col(data$hg38_fiveprime_breakpoint)
    data$hg38_threeprime_breakpoint <- trunc_break_col(data$hg38_threeprime_breakpoint)
    
    datatable(
      data,
      escape = FALSE, # NEW! escape = FALSE to show the hyperlink and clickable elipses
      selection = "single", # Allows one single row to be selected at a time
      rownames = FALSE,
      options = list(pageLength = 10, lengthMenu = c(10, 20, 50), scrollX = TRUE), # default will show 10 rows. User can pick to show 10, 20, or 50 rows. Horizontal scroll bar
      # click on expandable element (ellipsis) if there is one, then hide ellipsis and show other (hidden) breakpoints
      callback = JS("
      table.on('click', 'span.expand', function() {
        var $span = $(this);
        $span.hide();
        $span.siblings('span.full').show();
      });
    ")
    ) 
  })
  

  # Allow the user to download a tsv file with the fusions they filtered for (can download complete set if unfiltered) ####
  # Make a function now so that download_fusions_side and download_fusions_bottom output the same table
  make_fusion_download_handler <- function() {
    downloadHandler(
      filename = function() {
        paste("validatedFusions_", Sys.Date(), ".tsv", sep = "") # names file "validated fusions_year-month-date".tsv
      },
      content = function(file) {
        write.table(
          complete_filtered_df_to_download(),
          file,
          sep = "\t",
          row.names = FALSE,
          col.names = TRUE,
          quote = FALSE
        ) # will download the complete table based on what you have filtered (if no filtering, will download all entries). Not just entries displayed (10, 20, 50 per page), but all!!
      }
    )
  }
  
  # Run the same function for download_fusions_side download_fusions_bottom
  output$download_fusions_side <- make_fusion_download_handler()
  output$download_fusions_bottom <- make_fusion_download_handler()
  

  # Pull a selected row from "Table View" in prep for "Detail view" ####
  selected_row <- reactive({
    req(input$table_rows_selected) # ensure a row is selected
    complete_filtered_df_toDisplay() %>%
      slice(input$table_rows_selected) # use slice() to get row
  })
  
  # Top left card for "Detail view" - General Fusion Information ####
  output$topLeft_card <- renderUI({
    if (length(input$table_rows_selected) == 0) {
      return(
        tags$div(
          class = "alert alert-info",
          "Select a row from the 'Table View' tab to see details here."
        )
      )
    }

    row_data <- selected_row() %>%
      # Select relevant general information columns
      select(
        fusion_name_link,
        geneID_fusionName,
        cell_line_link,
        tissueType,
        discordant_read_support,
        breakpoint_support,
        program
      ) %>%
      # Change tissue types to be lower case, swap underscores for spaces, then uppercase first letter
      mutate(tissueType = str_replace_all(tolower(tissueType), "_", " ") %>%
        str_to_sentence()) %>%
      # Change supports to have capitalized first letter (lowercase rest)
      mutate(
        discordant_read_support = str_to_sentence(discordant_read_support),
        breakpoint_support = str_to_sentence(breakpoint_support)
      ) %>%
      # Update program display
      mutate(program = case_when(
        (program == "arriba") ~ "Arriba",
        (program == "star-fusion") ~ "STAR-Fusion",
        (program == "both") ~ "Arriba & STAR-Fusion" # if "one row" is select in side panel, then there will be instances of "both" in the DF
      )) %>%
      # Rename columns for better display in detail view
      rename(
        "Fusion" = "fusion_name_link",
        "Fusion (IDs)" = "geneID_fusionName",
        "Cell line" = "cell_line_link",
        "Tissue" = "tissueType",
        "Discordant read support" = "discordant_read_support",
        "Breakpoint support" = "breakpoint_support",
        "Fusion prediction program" = "program"
      )

    tags$div(
      class = "card",
      tags$div(
        style = "background-color: #834A88; color: white; padding: 10px; font-weight: bold;",
        "General Fusion Information"
      ),
      tags$div(
        class = "card-body",
        tags$table(
          class = "table table-borderless",
          tags$tbody(lapply(names(row_data), function(col_name) {
            value <- as.character(row_data[[col_name]])

            # If the value contains commas, split into new lines
            formatted_value <-
              if (grepl(",", value)) {
                # Split by commas and create line breaks
                formatted_value <-
                  paste0(strsplit(value, ",")[[1]], collapse = "<br>")
              } else {
                value
              }

            tags$tr(
              tags$td(tags$b(col_name)), # Column Name (Bold)
              tags$td(HTML(formatted_value)) # Displaying multiple lines with HTML <br> tags
            ) 
          }))
        )
      )
    )
  })
    
  # Top right card for "Detail view" - Additional Information ####
  output$topRight_card <- renderUI({
    if (length(input$table_rows_selected) == 0) {
      return(
        tags$div(
          class = "alert alert-info",
          "Select a row from the 'Table View' tab to see details here."
        )
      )
    }

    row_data <- selected_row() %>%
      # Select relevant additional information columns
      select(
        reading_frame_arriba,
        reading_frame_sf,
        program,
        fusionType_arriba,
        confidence_arriba,
        discordant_mates_arriba,
        junctionReadCount_sf,
        spanningFragCount_sf,
        ffpm_sf,
        promoterSwapEventCandidate,
        mitelmanFusion_match,
        mitelman_tissueTypes,
        tf_match,
        tf_disease,
        tcga_match,
        tcga_cancer_type
      ) %>%
      # Add a non-program specific call for reading frame
      mutate(
        reading_frame = case_when(
          (program == "arriba") ~ reading_frame_arriba,
          (program == "star-fusion") ~ reading_frame_sf,
          (program == "both" &
            reading_frame_arriba == reading_frame_sf) ~ reading_frame_arriba # if program is listed as both, then reading_frame_arriba ALWAYS is always the same as reading_frame_sf
        )
      ) %>%
      # Relocate the reading frame so it still shows up first
      relocate(reading_frame, .after = reading_frame_sf) %>%
      # De-select these cols that were only used to create reading_frame
      select(-reading_frame_arriba, -reading_frame_sf, -program) %>%
      # Update formatting to start with uppercase
      mutate(
        reading_frame = str_to_sentence(reading_frame),
        fusionType_arriba = str_to_sentence(fusionType_arriba),
        confidence_arriba = str_to_sentence(confidence_arriba),
        promoterSwapEventCandidate = str_to_sentence(promoterSwapEventCandidate),
        mitelmanFusion_match = str_to_sentence(mitelmanFusion_match),
        tf_match = str_to_sentence(tf_match),
        tcga_match = str_to_sentence(tcga_match),
      ) %>%
      
      # Fix mitelman tissue type formatting! Replace "_" and str_to_sentence tissue type separated by a comma
      mutate(mitelman_tissueTypes = strsplit(as.character(mitelman_tissueTypes), ",")) %>% # split by comma
      rowwise() %>%
      mutate(mitelman_tissueTypes = list(str_replace_all(mitelman_tissueTypes, "_", " "))) %>% # "_" to " "
      mutate(mitelman_tissueTypes = list(str_to_sentence(mitelman_tissueTypes))) %>% # uppercase first letter of tissue type
      mutate(mitelman_tissueTypes = paste(mitelman_tissueTypes, collapse = ",")) %>% # join again with a comma
      ungroup() %>%
      
      # Fix TF and TCGA disease formatting! Pull unique comma sep entries, replace "_" and str_to_sentence disease separated by a comma
      mutate(tf_disease = strsplit(as.character(tf_disease), ",")) %>% # split by comma
      rowwise() %>%
      mutate(tf_disease = list(unique(tf_disease))) %>% # pull unique only! comma sep in tsv for splitting purposes with the junctions
      mutate(tf_disease = list(str_replace_all(tf_disease, "_", " "))) %>% # "_" to " "
      mutate(tf_disease = list(str_to_sentence(tf_disease))) %>% # uppercase first letter of disease
      mutate(tf_disease = paste(tf_disease, collapse = ",")) %>% # join again with a comma
      ungroup() %>%
      
      mutate(tcga_cancer_type = strsplit(as.character(tcga_cancer_type), ",")) %>% # split by comma
      rowwise() %>%
      mutate(tcga_cancer_type = list(unique(tcga_cancer_type))) %>% # pull unique only! comma sep in tsv for splitting purposes with the junctions
      mutate(tcga_cancer_type = list(str_replace_all(tcga_cancer_type, "_", " "))) %>% # "_" to " "
      mutate(tcga_cancer_type = list(str_to_sentence(tcga_cancer_type))) %>% # uppercase first letter of cancer type
      mutate(tcga_cancer_type = paste(tcga_cancer_type, collapse = ",")) %>% # join again with a comma
      ungroup() %>%
      
      # String manipulation above affects NAs! Reset so that if it is "False" (remember we used str_to_sentence) that it resets tissue/disease to NA so that the upcoming NA filter catches it and doesn't display it
      mutate(mitelman_tissueTypes = if_else(mitelmanFusion_match == "False", NA, mitelman_tissueTypes),
             tf_disease = if_else(tf_match == "False", NA, tf_disease),
             tcga_cancer_type = if_else(tcga_match == "False", NA, tcga_cancer_type)) %>%
      
      rename(
        "Reading frame" = reading_frame,
        "Fusion type (Arriba)" = fusionType_arriba,
        "Confidence (Arriba)" = confidence_arriba,
        "Discordant mates (Arriba)" = discordant_mates_arriba,
        "Junction read count (STAR-Fusion)" = junctionReadCount_sf,
        "Spanning fragment count (STAR-Fusion)" = spanningFragCount_sf,
        "Fusion fragments per million (STAR-Fusion)" = ffpm_sf,
        "Promoter swapping candidate" = promoterSwapEventCandidate,
        "Mitelman fusion match" = mitelmanFusion_match,
        "Mitelman tissue types" = mitelman_tissueTypes,
        "TumorFusions match" = tf_match,
        "TumorFusions disease" = tf_disease,
        "TCGA match" = tcga_match,
        "TCGA cancer type" = tcga_cancer_type
      ) %>%
      # Only include columns that aren't empty
      select(where(~ all(!is.na(.)) & all(. != "")))

    tags$div(
      class = "card",
      tags$div(
        style = "background-color: #824C5F; color: white; padding: 10px; font-weight: bold;",
        "Additional Information"
      ),
      tags$div(
        class = "card-body",
        tags$table(
          class = "table table-borderless",
          tags$tbody(lapply(names(row_data), function(col_name) {
            value <- as.character(row_data[[col_name]])

            # If the value contains commas, split into new lines
            formatted_value <-
              if (grepl(",", value)) {
                # Split by commas and create line breaks
                formatted_value <-
                  paste0(strsplit(value, ",")[[1]], collapse = "<br>")
              } else {
                value
              }

            tags$tr(
              tags$td(tags$b(col_name)), # Column Name (Bold)
              tags$td(HTML(formatted_value)) # Displaying multiple lines with HTML <br> tags
            ) 
          }))
        )
      )
    )
  })
    
  # Bottom left card for "Detail view" - 5' partner Information ####
  output$bottomLeft_card <- renderUI({
    if (length(input$table_rows_selected) == 0) {
      return(
        tags$div(
          class = "alert alert-info",
          "Select a row from the 'Table View' tab to see details here."
        )
      )
    }

    row_data <- selected_row() %>%
      # Select 5' partner specific columns
      select(contains("fiveprime")) %>%

      # Change these to start with uppercase
      mutate(
        fiveprime_transcriptionType = str_to_sentence(fiveprime_transcriptionType),
        fiveprime_miRNA_host = str_to_sentence(fiveprime_miRNA_host),
        fiveprime_kinase = str_to_sentence(fiveprime_kinase),
        cgc_fiveprime_gene = str_to_sentence(cgc_fiveprime_gene)
      ) %>%
      # Change junction site to start with uppercase
      mutate(
        fiveprime_junctionSite_arriba = case_when(
          # Manually re-write instead of string to sentence because there are cases of all caps (e.g. "CDS")
          (fiveprime_junctionSite_arriba == "exon") ~ "Exon",
          (fiveprime_junctionSite_arriba == "exon/splice-site") ~ "Exon/splice-site",
          (fiveprime_junctionSite_arriba == "intron") ~ "Intron",
          TRUE ~ fiveprime_junctionSite_arriba # keeps the other scenarios as-is, e.g. "CDS", "5' UTR"
        )
      ) %>%
      # Update this column for clarity after removing the underscore
      mutate(fiveprime_inFrame_kinase = case_when(
        (fiveprime_inFrame_kinase == "FRAME_UNKNOWN") ~ "Unknown frame",
        TRUE ~ fiveprime_inFrame_kinase # keeps the other scenarios as-is, i.e. "TRUE" and "FALSE
      )) %>%
      # Get rid of underscore and start with uppercase
      mutate(
        fiveprime_gene_type = case_when(
          (fiveprime_gene_type == "protein_coding") ~ "Protein coding",
          (
            fiveprime_gene_type == "transcribed_unprocessed_pseudogene"
          ) ~ "Transcribed unprocessed pseudogene",
          (fiveprime_gene_type == "transcribed_processed_pseudogene") ~ "Transcribed processed pseudogene",
          (fiveprime_gene_type == "transcribed_unitary_pseudogene") ~ "Transcribed unitary pseudogene",
          (fiveprime_gene_type == "unprocessed_pseudogene") ~ "Unprocessed pseudogene",
          (fiveprime_gene_type == "artifact") ~ "Artifact",
          TRUE ~ fiveprime_gene_type # if lncRNA or TEC, will remain unchanged
        )
      ) %>%
      rename(
        "Gene type" = fiveprime_gene_type,
        "Gene coordinates" = hg38_fiveprime_gene_gtf_coordinates,
        "Junction coordinate" = hg38_fiveprime_junction_coordinate,
        "Genomic breakpoint" = hg38_fiveprime_breakpoint_sanityPassingPair,
        "Transcription type" = fiveprime_transcriptionType,
        "Split reads (Arriba)" = fiveprime_split_reads_arriba,
        "Coverage (Arriba)" = fiveprime_coverage_arriba,
        "Junction site (Arriba)" = fiveprime_junctionSite_arriba,
        "miRNA host" = fiveprime_miRNA_host,
        "miRNA hosted" = miRNA_fiveprime,
        "miRNA hosted - coordinates" = miRNA_fiveprime_coords,
        "Kinase partner" = fiveprime_kinase,
        "In-frame kinase" = fiveprime_inFrame_kinase,
        "Kinase" = fiveprime_kinase_name,
        "Kinase UniProt ID" = fiveprime_kinase_uniprotID,
        "Kinase group" = fiveprime_kinase_group,
        "Kinase family" = fiveprime_kinase_family,
        "Kinase sub-family" = fiveprime_kinase_subFamily,
        "CGC gene" = cgc_fiveprime_gene,
        "CGC description" = cgc_fiveprime_name,
        "CGC tier" = cgc_fiveprime_tier,
        "CGC hallmark" = cgc_fiveprime_hallmark,
        "CGC tumour types - somatic" = cgc_fiveprime_tumour_types_somatic,
        "CGC tumour types - germline" = cgc_fiveprime_tumour_types_germline,
        "CGC cancer syndrome" = cgc_fiveprime_cancer_syndrome,
        "CGC tissue type" = cgc_fiveprime_tissue_type,
        "CGC role in cancer" = cgc_fiveprime_role_in_cancer,
        "CGC mutation types" = cgc_fiveprime_mutation_types
      ) %>%
      # Only include columns that aren't empty
      select(where(~ all(!is.na(.)) & all(. != "")))

    tags$div(
      class = "card",
      tags$div(
        style = "background-color: #49895D; color: white; padding: 10px; font-weight: bold;",
        "5' Partner"
      ),
      tags$div(
        class = "card-body",
        tags$table(
          class = "table table-borderless",
          tags$tbody(lapply(names(row_data), function(col_name) {
            value <- as.character(row_data[[col_name]])

            # If the value contains commas, split into new lines
            formatted_value <-
              if (grepl(",", value)) {
                # Split by commas and create line breaks
                formatted_value <-
                  paste0(strsplit(value, ",")[[1]], collapse = "<br>")
              } else {
                value
              }

            tags$tr(
              tags$td(tags$b(col_name)), # Column Name (Bold)
              tags$td(HTML(formatted_value)) # Displaying multiple lines with HTML <br> tags
            ) 
          }))
        )
      )
    )
  })
  
  # Bottom right card for "Detail view" - 3' partner Information ####
  output$bottomRight_card <- renderUI({
    if (length(input$table_rows_selected) == 0) {
      return(
        tags$div(
          class = "alert alert-info",
          "Select a row from the 'Table View' tab to see details here."
        )
      )
    }

    row_data <- selected_row() %>%
      # Select 3' partner specific columns
      select(contains("threeprime")) %>%
      # Change these to start with uppercase
      mutate(
        threeprime_transcriptionType = str_to_sentence(threeprime_transcriptionType),
        threeprime_miRNA_host = str_to_sentence(threeprime_miRNA_host),
        threeprime_candidate_miRNA_host = str_to_sentence(threeprime_candidate_miRNA_host),
        threeprime_kinase = str_to_sentence(threeprime_kinase),
        cgc_threeprime_gene = str_to_sentence(cgc_threeprime_gene)
      ) %>%
      # Change junction site to start with uppercase
      mutate(
        threeprime_junctionSite_arriba = case_when(
          # manually re-write instead of string to sentence because there are cases of all caps (e.g. "CDS")
          (threeprime_junctionSite_arriba == "exon") ~ "Exon",
          (threeprime_junctionSite_arriba == "exon/splice-site") ~ "Exon/splice-site",
          (threeprime_junctionSite_arriba == "intron") ~ "Intron",
          TRUE ~ threeprime_junctionSite_arriba # keeps the other scenarios as-is
        )
      ) %>%
      # Update this column for clarity after removing the underscore
      mutate(threeprime_inFrame_kinase = case_when(
        (threeprime_inFrame_kinase == "FRAME_UNKNOWN") ~ "Unknown frame",
        TRUE ~ threeprime_inFrame_kinase # keeps the other scenarios as-is, i.e. "TRUE" and "FALSE
      )) %>%
      # Get rid of underscore and start with uppercase
      mutate(
        threeprime_gene_type = case_when(
          (threeprime_gene_type == "protein_coding") ~ "Protein coding",
          (
            threeprime_gene_type == "transcribed_unprocessed_pseudogene"
          ) ~ "Transcribed unprocessed pseudogene",
          (threeprime_gene_type == "transcribed_processed_pseudogene") ~ "Transcribed processed pseudogene",
          (threeprime_gene_type == "transcribed_unitary_pseudogene") ~ "Transcribed unitary pseudogene",
          (threeprime_gene_type == "unprocessed_pseudogene") ~ "Unprocessed pseudogene",
          (threeprime_gene_type == "processed_pseudogene") ~ "Processed pseudogene",
          TRUE ~ threeprime_gene_type # if lncRNA or TEC, will remain unchanged
        )
      ) %>%
      rename(
        "Gene type" = threeprime_gene_type,
        "Gene coordinates" = hg38_threeprime_gene_gtf_coordinates,
        "Junction coordinate" = hg38_threeprime_junction_coordinate,
        "Genomic breakpoint" = hg38_threeprime_breakpoint_sanityPassingPair,
        "Transcription type" = threeprime_transcriptionType,
        "Split reads (Arriba)" = threeprime_split_reads_arriba,
        "Coverage (Arriba)" = threeprime_coverage_arriba,
        "Junction site (Arriba)" = threeprime_junctionSite_arriba,
        "miRNA host" = threeprime_miRNA_host,
        "miRNA hosted" = miRNA_threeprime,
        "miRNA hosted - coordinates" = miRNA_threeprime_coords,
        "Candidate miRNA host" = threeprime_candidate_miRNA_host,
        "Downstream miRNA" = miRNA_downstream_threeprime,
        "Downstream miRNA - coordinates" = miRNA_downstream_threeprime_coords,
        "Kinase partner" = threeprime_kinase,
        "In-frame kinase" = threeprime_inFrame_kinase,
        "Kinase" = threeprime_kinase_name,
        "Kinase UniProt ID" = threeprime_kinase_uniprotID,
        "Kinase group" = threeprime_kinase_group,
        "Kinase family" = threeprime_kinase_family,
        "Kinase sub-family" = threeprime_kinase_subFamily,
        "CGC gene" = cgc_threeprime_gene,
        "CGC description" = cgc_threeprime_name,
        "CGC tier" = cgc_threeprime_tier,
        "CGC hallmark" = cgc_threeprime_hallmark,
        "CGC tumour types - somatic" = cgc_threeprime_tumour_types_somatic,
        "CGC tumour types - germline" = cgc_threeprime_tumour_types_germline,
        "CGC cancer syndrome" = cgc_threeprime_cancer_syndrome,
        "CGC tissue type" = cgc_threeprime_tissue_type,
        "CGC role in cancer" = cgc_threeprime_role_in_cancer,
        "CGC mutation types" = cgc_threeprime_mutation_types
      ) %>%
      # Only include columns that aren't empty
      select(where(~ all(!is.na(.)) & all(. != "")))

    tags$div(
      class = "card",
      tags$div(
        style = "background-color: #c85103; color: white; padding: 10px; font-weight: bold;",
        "3' Partner"
      ),
      tags$div(
        class = "card-body",
        tags$table(
          class = "table table-borderless",
          tags$tbody(lapply(names(row_data), function(col_name) {
            value <- as.character(row_data[[col_name]])

            # If the value contains commas, split into new lines
            formatted_value <-
              if (grepl(",", value)) {
                # Split by commas and create line breaks
                formatted_value <-
                  paste0(strsplit(value, ",")[[1]], collapse = "<br>")
              } else {
                value
              }

            tags$tr(
              tags$td(tags$b(col_name)), # Column Name (Bold)
              tags$td(HTML(formatted_value))
            ) # Displaying multiple lines with HTML <br> tags
          }))
        )
      )
    )
  })
  
  # "Cell Line Supplementary" tab section ####
  output$supp_blurb <- renderUI({
    if (input$tabs_ID == "Cell Line Supplementary") {
      # When "Cell Line Supplementary" is selected
      tagList(
        p(""),
        p("Here you will find more information about each CCLE cell line and its validated gene fusions."),
        p("To download the table, select 'Download supplementary (TSV)' at the bottom of your screen.")
      )
    } else {
      NULL # don't show this blurb for other tabs (i.e. no blurb for Detail View)
    }
  })
  
  output$supp_table <- renderDT({
    # Join the original supp_cellLine to replace the regular sampleID col with a cell line link. This is ONLY for the display version, NOT the download
    supp_cellLine_withLink <- supp_cellLine %>%
      left_join(cell_line_links) %>%
      relocate(cell_line_link, .after = sampleID) %>%
      select(-sampleID) %>%
      rename(sampleID = "cell_line_link")

    datatable(
      supp_cellLine_withLink,
      escape = FALSE,
      selection = "single",
      options = list(pageLength = 20, lengthMenu = c(10, 20, 50), scrollX = TRUE),
      rownames = FALSE
    ) # default will show 20 rows. User can pick to show 10, 20, or 50 rows. Allows one single row to be selected at a time. escape = FALSE for the hyperlink
  })
  
  # Allow the user to download a tsv file of the supplementary table ####
  output$download_supp <- downloadHandler(
    filename = function() {
      paste("cellLine_supplement.tsv") # names file "cellLine_supplement.tsv"
    },
    content = function(file) {
      write.table(
        supp_cellLine,
        file,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
      ) # will download the supplementary table as a tsv
    }
  )
  
  # "About" tab section ####
  output$col_blurb <- renderUI({
    if (input$tabs_ID == "About") {
      # when "About" is selected
      tagList(
        p(""),
        p("Here you will find descriptions of what each column contains."),
        p("For the most accurate information, please follow the link to the resource which the column originates from."),
      )
    } else {
      NULL # don't show this blurb for other tabs (i.e. no blurb for Detail View)
    }
  })

    
  output$more_infoTable <-
    renderUI({
      # Create an HTML table, needs to be html table so the url can be clickable!
      table_html <-
        "<table class='table table-striped table-bordered'>"
      table_html <- # add a header!
        paste0(
          table_html,
          "<thead><tr><th>Resource</th><th>Website</th></tr></thead><tbody>"
        )

      # Update the rows to html
      for (i in 1:nrow(col_moreInfo)) {
        table_html <- paste0(
          table_html,
          "<tr>",
          "<td>",
          col_moreInfo$Resource[i],
          "</td>",
          "<td><a href='",
          col_moreInfo$URL[i],
          "' target='_blank'>",
          col_moreInfo$URL[i],
          "</a></td>",
          "</tr>"
        )
      }

      table_html <-
        paste0(table_html, "</tbody></table>") # finish off the table

      HTML(table_html) # treat is as HTML, don't assume string
    })
  
  output$col_table <- renderDT({
    datatable(
      col_guide,
      selection = "single",
      options = list(pageLength = 5, lengthMenu = c(5, 10, 20, 50)),
      rownames = FALSE
    ) # default will show 10 rows. User can pick to show 5, 10, 20, or 50 rows. Allows one single row to be selected at a time
  })
  
}

shinyApp(ui, server)