

server <- function(input, output, session) {
  normcounts <- reactiveVal()
  dedata <- reactiveVal()
  meta <- reactiveVal()
  genes <- reactiveVal()
  species <- reactiveVal()
  species_choices <- reactiveVal()
  
  ##### Load Data
  observeEvent(input$loadData, {
    
    
    #first the normalized DE data
    normcounts <<- imap(list.files('/data/normcounts', full.names = T), read.csv)
    names(normcounts) <<- paste(str_split_i(list.files('/data/normcounts'), pattern = "_", i=1), 
                           str_split_i(list.files('/data/normcounts'), pattern = "_", i=2), sep='_')
    print(getwd())
    cnames <- c()
    specnames <- c()
    for(i in 1:length(normcounts)){
      cnames <- c(cnames, colnames(normcounts[[i]])[2:length(normcounts[[i]])])
      specnames <- c(specnames, rep(str_split_i(names(normcounts[i]),"_",1), length(normcounts[[i]]) -1 ))
    }
    
    #create meta from sample names
    meta <<- tibble(sample=cnames, species=specnames) %>% 
      mutate(individual = str_split_i(sample, "\\.|p", 1)) %>% 
      mutate(exp_group = case_when(
        str_sub(sample, -1, -1) == 'C' ~ "temperature",
        str_sub(sample, -1, -1) == 'H' ~ "hypoxia",
        str_sub(sample, -1, -1) == 'M' ~ "glucose",
        .default = NA
      )) %>% 
      mutate(treatment = case_when(
        str_sub(sample, -3, -1) == '5mM' ~ 2.5,
        str_sub(sample, -3, -1) == '8mM' ~ 8,
        str_sub(sample, -4, -1) == '30mM' ~ 30,
        str_sub(sample, -2, -1) == '6H' ~ 6,
        str_sub(sample, -3, -1) == '24H' ~ 24,
        str_sub(sample, -3, -1) == '32C' ~ 32,
        str_sub(sample, -3, -1) == '41C' ~ 41,
        .default = NA
      )) %>% 
      unique()
    
    #load DE data
    dedata <<- imap(list.files('/data/DESeq', full.names = T), read.csv)
    names(dedata) <<- apply( str_split_fixed(list.files('/data/DESeq'), pattern='_', n=4)[,1:3] , 1 , paste , collapse = "_" )
    
    
    #display DE Table
    output$summaryTable <- renderDT({
      detable <- imap_dfr(dedata, ~  data.frame(contrast= .y, nsig = filter(.x, padj<0.05) %>% nrow(), 
                                                nsig_up = filter(.x, padj<0.05, log2FoldChange>0) %>% nrow(),
                                                nsig_down = filter(.x, padj<0.05, log2FoldChange<0) %>% nrow()) ) %>% 
        separate_wider_delim(contrast, delim='_', names=c('species','exp_group','contrast'))
      DT::datatable(detable)
    })
    
    #create several UI outputs
    genes <<- normcounts %>% map(~.x$gene) %>% flatten_chr() %>% unique() 
    species <<- unique(meta$species)
    species_choices <<- setNames(species, str_to_title(species))
    
    output$geneSelect<- renderUI({
      selectizeInput(inputId='igene',
                  'Choose a Gene',
                  choices = genes,
                  multiple=F,
                  selected="TRIM32")
    })
    output$speciesSelectExp <- renderUI({
      selectizeInput("ispeciesexp", 
                   "Choose one or more species",
                   choices = species_choices, 
                   multiple=T, 
                   selected="human") 
    })
    output$speciesSelectDE <- renderUI({
      selectizeInput("ispeciesde", 
                     "Choose one or more species (Plot)",
                     choices = species_choices, 
                     multiple=T, 
                     selected="human") 
    })
    output$speciesSelectDETable <- renderUI({
      selectizeInput("ispeciesdetable", 
                     "Choose one species (Table)",
                     choices = species_choices, 
                     multiple=T, 
                     selected="human") 
    })
    
    
  })
  
  ##### Expression Plot
  observeEvent(input$refreshExpressionPlot, {
    #more or less need to output this directly still
    output$expressionPlot <- renderPlotly({ 
      req(input$ispeciesexp)
      req(input$iexpGroup)
      req(input$igene)
      
      ispecies <- input$ispeciesexp
      iexp_group <- input$iexpGroup
      igene <- input$igene

      isamps <- meta %>% filter(species %in% ispecies, exp_group==iexp_group) %>% pull(sample)
      datasets <- paste(ispecies, iexp_group, sep='_')
      expr_data <- normcounts[datasets] %>% 
        map(~ filter(.x, gene == igene)) %>%
        map(~ select(.x, any_of(isamps))) %>%
        compact() %>% 
        bind_cols() %>% 
        pivot_longer(
          cols = everything(),  
          names_to = "sample",    
          values_to = "expression" 
        )
      gdata <- left_join(expr_data, meta, by='sample') %>% 
        mutate(treatment = factor(treatment))
      
      #grouping variable color_axis
      x_axis <- input$iexpX
      if (x_axis=='individual') {
        #dot plot with treatment colors
        gg <- ggplot(gdata)+
          geom_point(aes(x=individual, y=expression, color=treatment, label=treatment))+
          facet_wrap(~species, scales='free_x', labeller = as_labeller(str_to_title) ) +
          theme_bw()+
          ylab("Normalized Expression (VST)")+
          xlab("Individual")+
          labs(color="Treatment")+
          ggtitle(paste0(igene, " ", str_to_title(iexp_group), " Expression in ", paste(str_to_title(ispecies), collapse=', ')))
        
      } else {
        #box plot with treatment colors
        gg <- ggplot(gdata, aes(x=treatment, y=expression, fill=treatment, label=individual))+
          geom_boxplot(outliers = F)+
          geom_jitter(height=0, width=0.1)+
          facet_wrap(~species, scales='free_x', labeller = as_labeller(str_to_title) ) +
          theme_bw()+
          ylab("Normalized Expression (VST)")+
          xlab("Treatment")+
          guides(fill="none")+
          ggtitle(paste0(igene, " ", str_to_title(iexp_group), " Expression in ", paste(str_to_title(ispecies), collapse=', ')))
        
      } 
       ggplotly(gg, tooltip = 'label')
      
    })
  })

  ##### DE Plot
  observeEvent(input$refreshDEPlot, {
    req(input$ideGroup)
    req(input$ispeciesde)
    req(input$deGuides)
    req(input$deColor)
    req(input$ispeciesdetable)
    
    iexp_group <- input$ideGroup
    ispecies <- input$ispeciesde
    deGuides <- input$deGuides
    deColor <- input$deColor
    
    
      if(iexp_group == "glucose") {
        icontrast <- input$glucoseContrast
      } else if(iexp_group == "hypoxia") {
        icontrast <- input$hypoxiaContrast
      } else if(iexp_group == "temperature") {
        icontrast <- input$tempContrast
      }
    
    
    
    datasets <- paste(ispecies, iexp_group, icontrast, sep='_')
    datasets <- datasets[datasets %in% names(dedata)] #ensure contrasts present in data
    
    #create a volcano plot
    lfccut <- 2
    pcut <- 0.05
    
    
    subset <- dedata[datasets] %>% 
      map(~ mutate(.x, sig=case_when(
        abs(log2FoldChange)>lfccut & padj < pcut ~ T,
        .default = F)) %>% 
          mutate(.x, expr_bin = cut(log10(.x$baseMean), 
                                breaks = quantile(log10(.x$baseMean), 
                                                  probs = seq(0, 1, 0.25), 
                                                  na.rm = TRUE),
                                labels = c("Low", "Medium-Low", "Medium-High", "High"),
                                include.lowest = TRUE))
        ) %>% 
      bind_rows(., .id='dataset')
    
    
    
    if (deColor=='expr'){icolor<-sym('expr_bin'); iclegend <- 'Expression Quartile'; fcolor <- scale_color_viridis_d()
    }else{icolor <- sym('sig');iclegend <- 'Significant'; fcolor <- scale_color_manual(values = c('black','red'))}
    
    #ouput plot
    output$dePlot <- renderPlotly({
        gg <- ggplot(subset)+
          geom_point(mapping=aes(x=log2FoldChange, y= - log10(padj), color=!!icolor, label=gene), alpha=0.7 )+
          fcolor+
          theme_bw()+
          facet_wrap(~dataset)+
          labs(color=iclegend, 
               title = paste0(str_to_title(iexp_group)," ", icontrast, " Volcano Plot in ", paste(str_to_title(ispecies), collapse=', ')))
        
        if (deGuides) { gg <- gg + 
            geom_hline(aes(yintercept=-log10(pcut)), color='red', alpha=0.8)+
            geom_vline(aes(xintercept=lfccut), color='grey', alpha=0.7)+
            geom_vline(aes(xintercept=-lfccut), color='grey', alpha=0.7)}
        
        ggplotly(gg, tooltip="label")
    
    })
    
    #output table
    output$deGenesTable <- renderDT({
      ispecies <- input$ispeciesdetable
      DT::datatable(dedata[[paste(ispecies, iexp_group, icontrast, sep='_')]])
    })
    
    
  })
  
}
