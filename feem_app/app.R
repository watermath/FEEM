#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load libraries
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinybusy)
library(shinyWidgets)

library(tidyverse)
library(ggplot2)
library(plotly)
library(cowplot)
library(ggrepel)
library(ggdendro)

library(lubridate)
library(psych)
library(corrplot)

# options(shiny.maxRequestSize=30*1024^2) 
# set.seed(1961)
theme_set(theme_bw())
appalette<-c("#000000","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
             "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

## FEEM functions ######
read_feem <- function(file) {
    check_file<-read_csv(file)
    if(dim(na.omit(check_file[1:2,]))[1]==0){
        samp_file<-read_csv(file)%>%
            slice(-c(1,2))%>%# remove first two rows from each FEEM after header
            pivot_longer(-1, names_to = "ex", values_to = "intensity") %>% 
            mutate_all(~ifelse(.x < 0, 0, .x) %>% as.numeric()) # set negative intensities to zero 
        
    }else{
    samp_file<-read_csv(file)%>%
        pivot_longer(-1, names_to = "ex", values_to = "intensity") %>% 
        mutate_all(~ifelse(.x < 0, 0, .x) %>% as.numeric()) # set negative intensities to zero
    }
    return(samp_file)
}

# Define FRI function
frinteg <- function(x, ex_min = 0, ex_max = 1e3, em_min = 0, em_max = 1e3) { 
    
    resol <- x %>% # average cell area
        summarize(
            ex_inc = unique(ex) %>% diff() %>% mean() %>% abs(),
            em_inc = unique(em) %>% diff() %>% mean() %>% abs(),
            fac = ex_inc * em_inc
        ) %>% 
        pull(fac)
    
    f_sum <- x %>% # total fluorescence intensity in region
        filter(
            ex > ex_min,
            ex <= ex_max,
            em > em_min,
            em <= em_max
        ) %>% 
        pull(intensity) %>% 
        sum()
    
    f_sum * resol # integrated area
}

scaleFUN <- function(x) sprintf("%.2f", x)

# sd_z<-function(x){(x - mean(x)) / sd(x)}
norm_z<-function(x){(x-min(x))/(max(x)-min(x))}

## FEEM panel
### User interface ###
ui <- dashboardPage(
    dashboardHeader(title = "FEEM Analysis"),
    
    # Sidebar content
    dashboardSidebar(
        
        sidebarMenu(
            menuItem("Plotting", tabName = "plotting", icon = icon("th")),
            menuItem("Other widgets", tabName = "other", icon = icon("dashboard"))
        ),
            
            fileInput("feemfiles", "Select zip folder of FEEM sample CSV files",multiple = FALSE),
                              
            fileInput("feem_metadata", 
                      span("Select FEEM metadata",
                           tags$a("(Required Format .csv)",
                                  href="#",
                                  onclick = "window.open('meta_input_format_example.png', 'newwindow'); return false;")
                      ),
                      accept = c(".csv", ".txt", ".tsv"),multiple = FALSE, placeholder = "No metadata here yet"
                              ),
                              
                          
        selectInput("disc_select",
                    label="Select categorical variable",
                    choices = "No metadata yet",
                    multiple = FALSE
                          ),    
        selectInput("disc_order",
                    label="Select categorical order",
                    choices = "No variable selected",
                    multiple = TRUE
        ), 
        
        selectInput("cts_select",
                    label="Select continuous variable",
                    choices = "No metadata yet",
                    multiple = FALSE
        ),  
                          # dateRangeInput("daterange",label="Select date range"),
                          
        numericInput("var_cutoff",label="Optional Ex/Em variance threshold",value = 0.1),
                          
        selectInput("feemsample",
                    label="Select individual FEEM sample",
                    choices = "No sample selected",
                    multiple = FALSE),
        
        actionButton("feemplot", label = "Plot sample FEEM",lib = "font-awesome", icon = icon("table-cells")),
                          
        actionButton("friplot", label="Plot FRI, peaks, indices",lib = "font-awesome", icon = icon("chart-column")),
        
        downloadBttn(
            outputId = "downloadRes",
            style = "bordered",
            color = "primary"
        )
        # downloadButton("downloadRes",label= "Download")
                          
                          
                          
    ),#end dashboardSidebar
             
    dashboardBody(
        
        tabItems(
            tabItem(tabName = "plotting",
                 # fluidRow(box(textOutput("test"))),
                 fluidRow(column(2),column(8, plotlyOutput("feem2d")),column(2)),
                 
                 fluidRow(column(12,plotlyOutput("fri_plot"))),
                 fluidRow(column(12,plotlyOutput("peak_plot"))),
                 fluidRow(column(12,plotlyOutput("ix_plot"))),
                 fluidRow(column(6,plotOutput("fri_corrs"))), 
                 fluidRow(column(12,plotlyOutput("feem_metabox")))
                 
             )#end tabItem
             
         ),#end endtabItems
         
tags$img(src = "cwrs_logo.png", style="height:200px; width:100%")
         
)#end dashboardBody










)#end dashboardPage


### Server logic ###
server<-function(input, output, session){
    cdata <- session$clientData
    rvs<-reactiveValues()
    #FEEM server logic
    
    # Loading data unzips and combines FEEM sample files into one data frame
    observeEvent(input$feemfiles,{
        withProgress(message="Processing FEEM files",detail="This may take a while depending on the number of sample files",{
            dir.name<-paste0(gsub(".zip","",input$feemfiles$name,ignore.case = T),"_raw_feem")
        unzip(input$feemfiles$datapath,overwrite = TRUE, junkpaths = TRUE
        ,exdir = dir.name
        )
            
        feem_raw <- list.files(dir.name,pattern = "*.csv", full.names = TRUE)%>%set_names()%>%map_dfr(read_feem, .id = "file")
        
        rvs$feem<-feem_raw%>%
            select(em = `Sample - Blank`, ex, intensity, file)%>%
            mutate(file=file%>%gsub(paste0(dir.name,"/"),"",.)%>%gsub(".csv","",.))
        
        # print(head(feem$file))
        # saveRDS(rvs$feem,"feem_combined.rds")
    })#end withProgress
    })
    
    #Load and select metadata
    observeEvent(input$feem_metadata, {
        feem.metadata<-read.csv(input$feem_metadata$datapath,header=T, check.names = FALSE)
        feem.metadata<-mutate(feem.metadata,file=as.character(file))%>%filter(file%in%rvs$feem$file)
        rvs$feem.metadata<-feem.metadata
        
        updateSelectInput(session, "disc_select", choices = names(feem.metadata))
        updateSelectInput(session, "cts_select", choices = names(feem.metadata))
        updateSelectInput(session, "feemsample", choices = unique(feem.metadata$sampleId))
        
        feem.long<-left_join(rvs$feem,select(feem.metadata,sampleId, file),by="file")
        rvs$feem.long<-feem.long
        
        #TEST
        # output$test<-renderDataTable(feem.long)
        
    #FEEM wide
        # feem.wide<-feem.long%>%
        #     select(-file)%>%
        #     pivot_wider(names_from = c("ex","em"), values_from = intensity, names_sep = "_")%>%
        #     column_to_rownames("sampleId")
        
        # feem.wide<-feem.wide%>%
        #     select(which(colSums(.)!=0)&which(!is.na(.)))%>%
        #     select(which(sapply(., var)>quantile(sapply(., var),input$var_cutoff)))
        
        # rvs$feem.wide<-feem.wide
        
    })#end observeEvent feem metadata
    
    #Select metadata discrete variable column
    observeEvent(input$disc_select,{
        rvs$disc_select<-input$disc_select
        if(!is.null(rvs$feem.metadata)){
        init.order<-unique(rvs$feem.metadata%>%pull(rvs$disc_select))
        updateSelectInput(session, "disc_order", choices = init.order)
        }else{
        updateSelectInput(session, "disc_order",choices = "No metadata yet", selected = character(0))
        }
    })
    
    observeEvent(input$disc_order,{
        rvs$disc_order<-input$disc_order
        #TEST
        # output$test<-renderText(rvs$disc_order)
    })
    
    #Select metadata continuous variable column
    observeEvent(input$cts_select,{
        rvs$cts_select<-input$cts_select
    })
    
    
    # Click plot of individual feem samples
    observeEvent(input$feemplot,{
        withProgress(message="Plotting FEEM sample heatmap",{
        single.feem<-rvs$feem%>%
            filter(file%in%as.character(input$feemsample))
        
        
            feem2d<-ggplot(single.feem) + 
                geom_raster(aes(ex, em, fill=intensity), interpolate = TRUE) +
                scale_fill_gradientn(colours=c("blue","orange","red"))+
                scale_x_continuous(name="Excitation Wavelength (nm)", breaks = seq(250,600,by=50))+
                scale_y_continuous(name="Emission Wavelength (nm)", breaks = seq(200,600,by=100))+
                ggtitle("")+
                theme_minimal()+
                theme(plot.title = element_text(size=14, face = "bold"),
                      axis.title = element_text(size = 14, face = "bold"),
                      strip.text = element_text(size = 12,color="black", face = "bold"),
                      strip.background =element_rect(fill="white"),
                      legend.title =element_text(size=12,face = "bold"),  
                      legend.text = element_text(size=12,face = "bold"), 
                      # legend.position = c(0.75,0.25), 
                      legend.key.size = unit(1, 'cm'), 
                      axis.text = element_text(size = 12, face = "bold"))
            output$feem2d<-renderPlotly({feem2d%>%ggplotly()})
    })#end withProgress  
    })
    
    
    # Click plot of feem FRI and sample metadata
    observeEvent((input$friplot),{
        withProgress(message="Plotting FRI, peaks, and indices vs metadata selected",detail="This may take a moment",{
            
        feem_peak_ix<-rvs$feem.long%>%
            group_by(sampleId)%>%
            summarize(
                      A=intensity[(abs(ex-260)==min(abs(ex-260)))&(abs(em-450)==min(abs(em-450)))],
                      B=intensity[(abs(ex-275)==min(abs(ex-275)))&(abs(em-310)==min(abs(em-310)))],
                      C=max(intensity[ex>=320&ex<=340&em>=410&em<=430]),
                      M=max(intensity[ex>=310&ex<=320&em>=380&em<=420]),
                      T=intensity[(abs(ex-275)==min(abs(ex-275)))&(abs(em-350)==min(abs(em-350)))],
                      bix = intensity[(abs(ex-310)==min(abs(ex-310)))&(abs(em-380)==min(abs(em-380)))] / # should be ex/em 310/380
                    max(intensity[(abs(ex-310)==min(abs(ex-310))) & em>=420 & em<=435]),
                    fi = intensity[(abs(ex-370)==min(abs(ex-370)))&(abs(em-470)==min(abs(em-470)))] /
                        intensity[(abs(ex-370)==min(abs(ex-370)))&(abs(em-520)==min(abs(em-520)))],
                      hix = sum(intensity[(abs(ex-255)==min(abs(ex-255))) & em >= 434 & em <= 480]) / sum(intensity[(abs(ex-255)==min(abs(ex-255))) & em >= 300 & em <= 344])
            )%>%
            # mutate(
            #     AT=A/T,
            #     CA=C/A,
            #     CM=C/M,
            #     CT=C/T,
            # )
            ungroup()
            
        feem_FRI <- rvs$feem.long%>%
                group_by(sampleId) %>%
                nest() %>%
                mutate(
                    total = map(data, frinteg),
                    protein1 = map(data, ~frinteg(.x, 200, 250, 200, 330)),
                    protein2 = map(data, ~frinteg(.x, 200, 250, 330, 380)),
                    fulvic = map(data, ~frinteg(.x, 200, 250, 380, 550)),  
                    microbial = map(data, ~frinteg(.x, 250, 340, 200, 380)),                 humic = map(data, ~frinteg(.x, 250, 400, 380, 550)),
                    cyano = map(data, ~frinteg(.x, 550, 600, 550, 620)),
                    pigment = map(data, ~frinteg(.x, 350, 450, 400, 621))
                )%>%
                unnest(c(total:pigment)) %>%
                select_if(~ !is.list(.x))%>%
                ungroup()
            
        feem.fri<-full_join(feem_FRI,feem_peak_ix,by="sampleId")#%>%
            # mutate(across(where(is.numeric), ~na_if(., Inf)), across(where(is.numeric), ~na_if(., -Inf)))%>% 
            # mutate(across(protein1:humic, ~./total))#divide by total
            # mutate(across(protein1:humic, norm_z))#normalize 
        fri.meta<-full_join(select(rvs$feem.metadata, sampleId,rvs$disc_select,rvs$cts_select),select(feem.fri,-total), by="sampleId")
        
        rvs$fri.meta<-fri.meta
        
        #Conventional regions + cyano, pigment FRI barplot 
        fri.plot<-pivot_longer(fri.meta,cols=protein1:humic,names_to = "Region", values_to = "Volume")
        
        fri_plot<-ggplot(fri.plot,aes(x=factor(.data[[rvs$disc_select]], levels = rvs$disc_order)
            ,y=Volume, colour=Region
        ))+
            geom_boxplot()+
                # geom_point(size=3, position = position_jitter(height = 0))+
                # geom_bar(stat="identity",position=position_dodge())+
                # geom_text(size = 3, position = position_stack(vjust = 0.5))+
                scale_colour_manual(values=appalette)+
                # scale_fill_manual(values=appalette)+
                theme(plot.title = element_text(size=16, face = "bold"),
                      axis.title = element_text(size = 16, face = "bold"),
                      legend.title = element_text(size=16,face = "bold"),
                      legend.text = element_text(size=12,face = "bold"),
                      legend.position="bottom",
                      # axis.title.x=element_blank(),
                      # axis.text.x=element_blank(),
                      axis.text.x = element_text(angle=45, hjust=1),
                      axis.text = element_text(size = 16, face = "bold"))+
                # scale_y_continuous(labels = scales::percent)+
                ylab("FRI Volume")+xlab(rvs$disc_select)
        
        if(is.Date(rvs$disc_select)){
            fri_plot<-fri_plot+scale_x_date(
                # date_breaks = "1 months", date_labels = "%b-%y"
            )}
        
        output$fri_plot<-renderPlotly({fri_plot%>%ggplotly()%>%layout(boxmode = "group")})
        
        #Peak barplot
        peak.plot<-pivot_longer(fri.meta,cols=A:T,names_to = "Peak", values_to = "Intensity")
        peak_plot<-ggplot(peak.plot,aes(x=factor(.data[[rvs$disc_select]], levels = rvs$disc_order),y=Intensity, colour=Peak
        ))+
            geom_boxplot()+
            # geom_point(size=3, position = position_jitter(height = 0))+
            # geom_bar(stat="identity",position=position_dodge())+
            # geom_text(size = 3, position = position_stack(vjust = 0.5))+
            scale_colour_manual(values=appalette)+
            # scale_fill_manual(values=appalette)+
            theme(plot.title = element_text(size=16, face = "bold"),
                  axis.title = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size=16,face = "bold"),
                  legend.text = element_text(size=12,face = "bold"),
                  legend.position="bottom",
                  # axis.title.x=element_blank(),
                  # axis.text.x=element_blank(),
                  axis.text.x = element_text(angle=45, hjust=1),
                  axis.text = element_text(size = 16, face = "bold"))+
            # scale_y_continuous(labels = scales::percent)+
            ylab("Intensity")+xlab(rvs$disc_select)
        
        if(is.Date(rvs$disc_select)){
            peak_plot<-peak_plot+scale_x_date(
                # date_breaks = "1 months", date_labels = "%b-%y"
            )}
        
        output$peak_plot<-renderPlotly({peak_plot%>%ggplotly()%>%layout(boxmode = "group")})
        
        #Index barplot
        ix.plot<-pivot_longer(fri.meta,cols=bix:fi,names_to = "Index", values_to = "Intensity")
        ix_plot<-ggplot(ix.plot,aes(x=factor(.data[[rvs$disc_select]], levels = rvs$disc_order),y=Intensity, colour=Index
        ))+
            geom_boxplot()+
            # geom_point(size=3, position = position_jitter(height = 0))+
            # geom_bar(stat="identity",position=position_dodge())+
            # geom_text(size = 3, position = position_stack(vjust = 0.5))+
            scale_colour_manual(values=appalette)+
            # scale_fill_manual(values=appalette)+
            theme(plot.title = element_text(size=16, face = "bold"),
                  axis.title = element_text(size = 16, face = "bold"),
                  legend.title = element_text(size=16,face = "bold"),
                  legend.text = element_text(size=12,face = "bold"),
                  legend.position="bottom",
                  # axis.title.x=element_blank(),
                  # axis.text.x=element_blank(),
                  axis.text.x = element_text(angle=45, hjust=1),
                  axis.text = element_text(size = 16, face = "bold"))+
            # scale_y_continuous(labels = scales::percent)+
            ylab("Intensity ratio")+xlab(rvs$disc_select)
        
        if(is.Date(rvs$disc_select)){
            ix_plot<-ix_plot+scale_x_date(
                # date_breaks = "1 months", date_labels = "%b-%y"
            )}
        
        output$ix_plot<-renderPlotly({ix_plot%>%ggplotly()%>%layout(boxmode = "group")})
        
        #metadata covariate plot(s)
        feem_metabox<-ggplot(fri.meta,aes(x=factor(.data[[rvs$disc_select]], levels = rvs$disc_order),y=.data[[rvs$cts_select]]))+
                geom_boxplot()+ #geom_jitter(width = 0.2)+
                theme(plot.title = element_text(size=16, face = "bold"),
                      axis.title = element_text(size = 16, face = "bold"),
                      legend.title = element_text(size=16,face = "bold"),
                      legend.text = element_text(size=12,face = "bold"),
                      axis.text.x = element_text(angle=45, hjust=1),
                      axis.text = element_text(size = 16, face = "bold")
                )+ylab(rvs$cts_select)+xlab(rvs$disc_select)
                # scale_y_continuous(labels=scaleFUN)+xlab(rvs$disc_select)
        output$feem_metabox<-renderPlotly({feem_metabox%>%ggplotly()})#end metadata covariate boxplot
        
        
        # output$fribar<-renderPlot({plot_grid(fri.barplot,newfri.barplot,cov.plot,ncol = 1) })
        
        fri.corrs<-fri.meta%>%select(where(is.numeric))
        corrs<-corr.test(fri.corrs, method = "spearman")
        output$fri_corrs<-renderPlot({
            corrplot(corrs$r, diag = FALSE, type = "upper", sig.level = 0.01, method = "square", tl.srt=45, 
                     p.mat = corrs$p, 
                     insig = "blank", na.label = ".", tl.col = "black",col=colorRampPalette(c("blue","white","red"))(200))
            
        })
        
    })#end withProgress
    })
    
    output$downloadRes <- downloadHandler(
        filename = function() {
            paste("fri_peaks_indices_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            write.csv(rvs$fri.meta, file)
        }
    )
    
    
    
} #end server

shinyApp(ui=ui,server=server)
