######################################
# GWAS Mapping Precision (gwasMP)
# @Author zhili zheng <zhilizheng@uq.edu.au>;
# @Author Yang Wu <y.wu2@uq.edu.au>;
# @Author Jian Yang <jian.yang@uq.edu.au>;
# @Bugs report to Jian Yang <jian.yang@uq.edu.au>
# @Updated Feb.26 2017
# @Licence GPL v3
######################################

shinyUI(fluidPage(
  titlePanel("GWAS Mapping Precision"), 
  withMathJax(),
  sidebarLayout(
    sidebarPanel(
      helpText("This tool is designed to quantify the physical distance, linkage disequilibrium (LD) and difference in minor allele frequency (MAF) between the top associated SNPs identified from GWAS and the underlying causal variants. The results are from simulations based on whole-genome sequencing data (Wu et al. 2017). Shown is proportion of the GWAS top SNPs of which the corresponding causal variants are within a certain distance, LD and/or MAF difference threshold as specified in the options below."),
    
      selectInput("maf", 
        label = "Please select the MAF range of the variants you are interested",
        choices = c("Common (MAF > 0.01)", "Rare (0.0004 < MAF < 0.01)", "0.01 < MAF < 0.05", "0.05 < MAF < 0.10", 
                    "0.10 < MAF < 0.20",  "0.20 < MAF < 0.30", "0.30 < MAF < 0.40","0.40 < MAF < 0.50"),
        selected = "Common (MAF > 0.01)"
      ),

      hr(),

      p("Please specify the thresholds"),

      numericInput("pd", 
        label = "Physical distance (Kb) between causal variants and GWAS hits",
        min = 5, max = 1000, value = 100),
      
      numericInput("ld", 
        label = "LD r\\(^2\\) between causal variants and GWAS hits",
        min = 0.1, max = 1, value = 0.8),
      
      numericInput("maf_diff", label = "MAF difference between causal variants and GWAS hits",
        min = 0, max = 0.5, value = 0.05),
      
      HTML("<p><strong>Citation:</strong> Wu et al. (2017) Quantifying the mapping precision of genome-wide association studies using whole-genome sequencing data. Genome Biology, 18(1): 86. </p>
      <p><strong>Credits:</strong> <a href='mailto:zhilizheng@uq.edu.au'>Zhili Zheng</a>, <a href='mailto:y.wu2@uq.edu.au'>Yang Wu</a> and <a href=' http://researchers.uq.edu.au/researcher/2713'>Jian Yang</a> at the Program in <a href='http://cnsgenomics.com/'>Complex Trait Genomics</a>, The University of Queensland</p>
      <p><strong>Troubleshooting:</strong> <a href='mailto:jian.yang@uq.edu.au'>Jian Yang</a></p>")
    ),
  
    mainPanel(fluidPage(
      plotOutput("plot", height=290),
	    hr(),

      fluidRow(
        column(12, 
          tableOutput("prop")
        ),
      
        column(12,
          tableOutput("both_prop")
        )
        
      ),
      fluidRow(
        column(12,
          HTML("<div style='color:gray; font-size:85%'><p><strong>Note:</strong><br/>
          WGS: GWAS using whole genome sequencing data.<br/>1KPG3: GWAS using SNP array data imputed to 1KGP Phase 3.<br/>1KPG1: GWAS using SNP array data imputed to the 1000 Genome Project (1KGP) Phase 1.<br/>HAPMAP2: GWAS using SNP array data imputed to the HapMap Project Phase 2.</p></div>")
        )
      )
    ))
  )
))
