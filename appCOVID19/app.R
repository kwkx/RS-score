library(shiny)
library(shinydashboard)
library(rintrojs)
library(shinyhelper)
library(magrittr)
library(nnet)
library(ggplot2)
library(ggpubr)
#load("model.Rdata")
header <- dashboardHeader(title = "RS score")

sidebar <- dashboardSidebar(
  width=300,
  introjsUI(),
  introBox(
    fluidRow(
      column(12, 
             align = "middle",
             flowLayout(actionButton("button1", "Click for instructions")
             )
      )
      
    )


  ),
  fluidRow(
    column(12, 
           align = "middle",
           flowLayout(helpText(h4("Please upload expression - data below."))%>%
                        helper(icon = "question-circle",
                               colour = "white",
                               type = "markdown",
                               content = "mymarkdown")
           )
    )
    
  ),
  introBox(
    data.step = 1,
    data.intro = "Please input the 33 gene expression data of your patients.",
    fileInput("file1", "Upload CSV File",
              multiple = TRUE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"))
        
      
    
  ),
  introBox(
    data.step = 2,
    data.intro = "You can click here to see our example data. We recommend you to download our example data and then upload it as your first try!",
    fluidRow(
      column(12, 
             align = "middle",
             flowLayout(downloadButton('downloadData', 'Click to download example data', class = "butt"),
                        tags$head(tags$style(".butt{background-color:#000000;}"))
             )
      )

    )
    #downloadButton('downloadData', 'Click to download sampledata', class = "butt"),
    #tags$style(type='text/css', "#run_report { width:50%; margin-left: 100px;}"),
    #tags$head(tags$style(".butt{background-color:#000000;}"))
  )
  
  
  
)

body <- dashboardBody(
      HTML('<meta name="viewport" content="width=1024">'),
  
   
      width = 6,
      hr(),
      hr(),
      introBox(
        data.step = 3,
        data.intro = "Please input a valid csv file to see the RS score of this COVID-19 patient.",
        h3("When you input a valid expression data (CSV format), you will get the RS score of this COVID-19 patient here:"),
        textOutput("riskscore"),
        tags$head(tags$style("#riskscore{color: red;
                                 font-size: 20px;
                                 font-style: bold;
                                 }"
        )
        )
        
        
        
        
      ),


      
      introBox(
        data.step = 4,
        data.intro = "You can get your RS score level here which can be easily compared with other type severity score in different datasets.",
        h3("The raw figure of different severe score types among 5 datasets:"),

        plotOutput(outputId = "rawsevereplot",width = 750,height = 800),
        downloadButton("downloadData1", "Download"),
        h3("When you input a valid expression data (CSV format), the RS score level will be shown as red line below:"),
        plotOutput(outputId = "severeplot",width = 750,height = 800),
        downloadButton("downloadData2", "Download")

        
        
        
      ),
      introBox(
        data.step = 5,
        data.intro = "Our contact is here.",
        hr(),
        h4("If you have any problems with using the web-app, or suggestions for improvement, please contact Weikaixin Kong at 1510307407@pku.edu.cn"),
        hr(),
        h4("Department of Molecular and Cellular Pharmacology"),
        h4("School of Pharmaceutical Sciences"),
        h4("Peking University Health Science Center"),
        h4("100191 Beijing"),
        h4("China")
      )
      
    

)

ui <- dashboardPage(header, sidebar, body)

server <- function(input, output,session) {
  
  df <- reactive({
    df <- read.csv("severedara.csv",
                   header = T,
                   sep = ",",row.names = 1)
    df
  })
  
  output$downloadData1 <- downloadHandler(
    filename = function(){
      paste0("rawfigure", '.pdf')
    },
    content = function(file){
      
      bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#0D595B","#570D5B","#0D5B10")
      bioCol=bioCol[1:length(table(df()$group))]
      my_comparisons=list()
      my_comparisons[[1]]=c("Short hospital-free time (GSE177477)","Long hospital-free time (GSE177477)")
      my_comparisons[[2]]=c("Symptomatic (GSE157103)","Asymptomatic (GSE157103)")
      
      my_comparisons[[3]]=c("Critical (GSE172114)","Non-critical (GSE172114)")
      
      my_comparisons[[4]]=c("ICU (GSE152418)","Severe (GSE152418)")
      my_comparisons[[5]]=c("ICU (GSE152418)","Moderate (GSE152418)")
      my_comparisons[[6]]=c("Severe (GSE152418)","Moderate (GSE152418)")
      
      my_comparisons[[7]]=c("severity: 0 (GSE155454)","severity: 1 (GSE155454)")
      my_comparisons[[8]]=c("severity: 0 (GSE155454)","severity: 2 (GSE155454)")
      my_comparisons[[9]]=c("severity: 1 (GSE155454)","severity: 2 (GSE155454)")
      
      
      boxplot2=ggboxplot(df(), x="group", y="Score", fill="group",
                         xlab="Different datasets",
                         legend.title="Datasets",
                         ylab="RS score",
                         palette=bioCol
      )+ theme(legend.position = "none")+  geom_jitter(width =0.2)+
        stat_compare_means(comparisons=my_comparisons,method= "t.test")+
        rotate_x_text(90)

      pdf(file,width=12,height=12)
      print(boxplot2)
      dev.off()
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function(){
      paste0("figure-with-patient-severity", '.pdf')
    },
    content = function(file){
      bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#0D595B","#570D5B","#0D5B10")
      bioCol=bioCol[1:length(table(df()$group))]
      my_comparisons=list()
      my_comparisons[[1]]=c("Short hospital-free time (GSE177477)","Long hospital-free time (GSE177477)")
      my_comparisons[[2]]=c("Symptomatic (GSE157103)","Asymptomatic (GSE157103)")
      
      my_comparisons[[3]]=c("Critical (GSE172114)","Non-critical (GSE172114)")
      
      my_comparisons[[4]]=c("ICU (GSE152418)","Severe (GSE152418)")
      my_comparisons[[5]]=c("ICU (GSE152418)","Moderate (GSE152418)")
      my_comparisons[[6]]=c("Severe (GSE152418)","Moderate (GSE152418)")
      
      my_comparisons[[7]]=c("severity: 0 (GSE155454)","severity: 1 (GSE155454)")
      my_comparisons[[8]]=c("severity: 0 (GSE155454)","severity: 2 (GSE155454)")
      my_comparisons[[9]]=c("severity: 1 (GSE155454)","severity: 2 (GSE155454)")
      
      
      h=as.numeric(datariskscore())
      boxplot2=ggboxplot(df(), x="group", y="Score", fill="group",
                         xlab="Different datasets",
                         legend.title="Datasets",
                         ylab="RS score",
                         palette=bioCol
      )+ theme(legend.position = "none")+  geom_jitter(width =0.2)+
        stat_compare_means(comparisons=my_comparisons,method= "t.test")+
        geom_hline(aes(yintercept=h),color="red",size=1) +
        scale_y_continuous(breaks = sort(c(seq(min(df()$Score), max(df()$Score), length.out=5), h)))+rotate_x_text(90)
      
      
      pdf(file,width=12,height=12)
      print(boxplot2)
      dev.off()
    }
  )
  
  
  
  
  output$rawsevereplot <- renderPlot({
    
    
    bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#0D595B","#570D5B","#0D5B10")
    bioCol=bioCol[1:length(table(df()$group))]
    my_comparisons=list()
    my_comparisons[[1]]=c("Short hospital-free time (GSE177477)","Long hospital-free time (GSE177477)")
    my_comparisons[[2]]=c("Symptomatic (GSE157103)","Asymptomatic (GSE157103)")
    
    my_comparisons[[3]]=c("Critical (GSE172114)","Non-critical (GSE172114)")
    
    my_comparisons[[4]]=c("ICU (GSE152418)","Severe (GSE152418)")
    my_comparisons[[5]]=c("ICU (GSE152418)","Moderate (GSE152418)")
    my_comparisons[[6]]=c("Severe (GSE152418)","Moderate (GSE152418)")
    
    my_comparisons[[7]]=c("severity: 0 (GSE155454)","severity: 1 (GSE155454)")
    my_comparisons[[8]]=c("severity: 0 (GSE155454)","severity: 2 (GSE155454)")
    my_comparisons[[9]]=c("severity: 1 (GSE155454)","severity: 2 (GSE155454)")
    
    
    boxplot2=ggboxplot(df(), x="group", y="Score", fill="group",
                       xlab="Different datasets",
                       legend.title="Datasets",
                       ylab="RS score",
                       palette=bioCol
    )+ theme(legend.position = "none")+  geom_jitter(width =0.2)+
      stat_compare_means(comparisons=my_comparisons,method= "t.test")+
      rotate_x_text(90)
    
    boxplot2
    
    
    
    
    
  })
  
  output$severeplot <- renderPlot({
    

        bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#0D595B","#570D5B","#0D5B10")
        bioCol=bioCol[1:length(table(df()$group))]
        my_comparisons=list()
        my_comparisons[[1]]=c("Short hospital-free time (GSE177477)","Long hospital-free time (GSE177477)")
        my_comparisons[[2]]=c("Symptomatic (GSE157103)","Asymptomatic (GSE157103)")
        
        my_comparisons[[3]]=c("Critical (GSE172114)","Non-critical (GSE172114)")
        
        my_comparisons[[4]]=c("ICU (GSE152418)","Severe (GSE152418)")
        my_comparisons[[5]]=c("ICU (GSE152418)","Moderate (GSE152418)")
        my_comparisons[[6]]=c("Severe (GSE152418)","Moderate (GSE152418)")
        
        my_comparisons[[7]]=c("severity: 0 (GSE155454)","severity: 1 (GSE155454)")
        my_comparisons[[8]]=c("severity: 0 (GSE155454)","severity: 2 (GSE155454)")
        my_comparisons[[9]]=c("severity: 1 (GSE155454)","severity: 2 (GSE155454)")
        
        
        h=as.numeric(datariskscore())
        boxplot2=ggboxplot(df(), x="group", y="Score", fill="group",
                          xlab="Different datasets",
                          legend.title="Datasets",
                          ylab="RS score",
                          palette=bioCol
        )+ theme(legend.position = "none")+  geom_jitter(width =0.2)+
          stat_compare_means(comparisons=my_comparisons,method= "t.test")+
        geom_hline(aes(yintercept=h),color="red",size=1) +
         scale_y_continuous(breaks = sort(c(seq(min(df()$Score), max(df()$Score), length.out=5), h)))+rotate_x_text(90)

        boxplot2
     


    
    
  })
  
  
  
  
  
  observe_helpers(help_dir = "helpfiles")
  
  
  datasetInput <- reactive({
    df <- read.csv("example.csv",
                   header = T,
                   sep = ",")
    df
  })
  
  
  coxweight <- reactive({
    df3 <- read.csv("coxweight.csv",
                    header = T,row.names=1,
                    sep = ",")
    df3
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("example data", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = FALSE,col.names = F)
    }
  )
  
  
  
  

  
  datariskscore <- reactive({
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = T,
                   sep = ",",
                   row.names = 1)
    
    scorevector=vector(length=36)
    scorevector[1]=ifelse(df["ADAMTS3",1]>df["FAM110C",1], 1, 0)
    
    scorevector[2]=ifelse(df["ADAMTS3",1]>df["NTN4",1], 1, 0)
    
    scorevector[3]=ifelse(df["AMOT",1]>df["HLF",1], 1, 0)
    
    scorevector[4]=ifelse(df["AMOT",1]>df["NTN4",1], 1, 0)
    
    scorevector[5]=ifelse(df["B3GALT2",1]>df["CACNA2D3",1], 1, 0)
    
    scorevector[6]=ifelse(df["B3GALT2",1]>df["DUSP4",1], 1, 0)
    
    scorevector[7]=ifelse(df["B3GALT2",1]>df["NTN4",1], 1, 0)
    
    scorevector[8]=ifelse(df["CACNA2D3",1]>df["CYSLTR2",1], 1, 0)
    scorevector[9]=ifelse(df["CACNA2D3",1]>df["ZNF683",1], 1, 0)
    
    scorevector[10]=ifelse(df["CYSLTR2",1]>df["IL5RA",1], 1, 0)
    scorevector[11]=ifelse(df["CYSLTR2",1]>df["SIGLEC8",1], 1, 0)
    
    scorevector[12]=ifelse(df["ELOVL4",1]>df["FGFR2",1], 1, 0)
    scorevector[13]=ifelse(df["ELOVL4",1]>df["SLC4A10",1], 1, 0)
    
    scorevector[14]=ifelse(df["FAM110C",1]>df["WNT7A",1], 1, 0)
    
    scorevector[15]=ifelse(df["FGFR2",1]>df["HLF",1], 1, 0)
    
    scorevector[16]=ifelse(df["GPR174",1]>df["IDO1",1], 1, 0)
    scorevector[17]=ifelse(df["GPR174",1]>df["PCOLCE2",1], 1, 0)
    
    scorevector[18]=ifelse(df["GRAPL",1]>df["MTUS1",1], 1, 0)
    
    scorevector[19]=ifelse(df["HLF",1]>df["NTN4",1], 1, 0)
    scorevector[20]=ifelse(df["HLF",1]>df["OLFM1",1], 1, 0)
    
    scorevector[21]=ifelse(df["IL5RA",1]>df["PCOLCE2",1], 1, 0)
    scorevector[22]=ifelse(df["IL5RA",1]>df["PID1",1], 1, 0)
    scorevector[23]=ifelse(df["IL5RA",1]>df["TBC1D4",1], 1, 0)
    
    scorevector[24]=ifelse(df["ISM1",1]>df["NTN4",1], 1, 0)
    scorevector[25]=ifelse(df["ISM1",1]>df["SIGLEC8",1], 1, 0)
    
    scorevector[26]=ifelse(df["KIAA1671",1]>df["WNT7A",1], 1, 0)
    
    scorevector[27]=ifelse(df["KLRC2",1]>df["PRSS33",1], 1, 0)
    
    scorevector[28]=ifelse(df["KLRC3",1]>df["PCOLCE2",1], 1, 0)
    
    scorevector[29]=ifelse(df["LRRN3",1]>df["NEFL",1], 1, 0)
    scorevector[30]=ifelse(df["LRRN3",1]>df["PLXDC1",1], 1, 0)
    
    scorevector[31]=ifelse(df["MS4A2",1]>df["NEFL",1], 1, 0)
    
    scorevector[32]=ifelse(df["NEFL",1]>df["SIGLEC8",1], 1, 0)
    scorevector[33]=ifelse(df["NEFL",1]>df["SLC4A10",1], 1, 0)
    
    scorevector[34]=ifelse(df["PCOLCE2",1]>df["PID1",1], 1, 0)
    scorevector[35]=ifelse(df["PCOLCE2",1]>df["SIGLEC8",1], 1, 0)
    scorevector[36]=ifelse(df["PLXDC1",1]>df["ZNF683",1], 1, 0)

  
    scorefinal=scorevector[1]*coxweight()[1,1]+scorevector[2]*coxweight()[2,1]+scorevector[3]*coxweight()[3,1]+
                     scorevector[4]*coxweight()[4,1]+scorevector[5]*coxweight()[5,1]+scorevector[6]*coxweight()[6,1]+
                     scorevector[7]*coxweight()[7,1]+scorevector[8]*coxweight()[8,1]+scorevector[9]*coxweight()[9,1]+
                     scorevector[10]*coxweight()[10,1]+scorevector[11]*coxweight()[11,1]+scorevector[12]*coxweight()[12,1]+
                     scorevector[13]*coxweight()[13,1]+scorevector[14]*coxweight()[14,1]+scorevector[15]*coxweight()[15,1]+
                     scorevector[16]*coxweight()[16,1]+scorevector[17]*coxweight()[17,1]+scorevector[18]*coxweight()[18,1]+
                     scorevector[19]*coxweight()[19,1]+scorevector[20]*coxweight()[20,1]+scorevector[21]*coxweight()[21,1]+
                     scorevector[22]*coxweight()[22,1]+scorevector[23]*coxweight()[23,1]+scorevector[24]*coxweight()[24,1]+
                     scorevector[25]*coxweight()[25,1]+scorevector[26]*coxweight()[26,1]+scorevector[27]*coxweight()[27,1]+
                     scorevector[28]*coxweight()[28,1]+scorevector[29]*coxweight()[29,1]+scorevector[30]*coxweight()[30,1]+
                     scorevector[31]*coxweight()[31,1]+scorevector[32]*coxweight()[32,1]+scorevector[33]*coxweight()[33,1]+
                     scorevector[34]*coxweight()[34,1]+scorevector[35]*coxweight()[35,1]+scorevector[36]*coxweight()[36,1]-1.784271167
    scorefinal=round(scorefinal,digits = 2)
    as.numeric(scorefinal)
  })
  
  output$riskscore=renderText(
    {

        as.character(paste0("RS score = ",datariskscore()))
      
    }
  )
  
  
  
  
  observeEvent(input$button1,
               introjs(session, options = list("nextLabel"="Next",
                                               "prevLabel"="Prev",
                                               "skipLabel"="End tour")
               )
  )
  
  
  
}

shinyApp(ui, server)
