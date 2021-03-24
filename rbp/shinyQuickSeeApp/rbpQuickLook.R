library(R6)
library(shiny)
library(DataTableWidget)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
buttonStyle <- "margin: 5px; margin-right: 0px; font-size: 14px;"

DataTableDemoApp = R6Class("app",

    #--------------------------------------------------------------------------------
    private = list(dtw = NULL,
                   randomTextButton = NULL,
                   tbl = NULL,
                   msgBox = NULL),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(){
            printf("initializing demo")
            #private$msgBox =  msgBoxWidget$new(id="box1", title="table selection", boxWidth=600)
            private$tbl = get(load("tbl.gata2.RData"))[, c(1:8, 10:12, 9)]
            private$tbl$score <- round(private$tbl$score, digits=2)
            private$dtw = dataTableWidget$new(id="tbl.1", private$tbl,
                                              width="98%", height="1000px",
                                              pageLength=50,
                                              border="1px blue solid; border-radius: 10px;")
                                              #columnWidths=c(50,60,70,80,90,100,110,120))
            },

        #------------------------------------------------------------
        ui = function(){
           fluidPage(
               private$dtw$ui()
              )},

        #------------------------------------------------------------
        server = function(input, output, session){

            printf("entering dataTableWidgetDemo::server")
            #private$msgBox$server(input, output, session)
            private$dtw$server(input, output, session)

            #observe({
            #   mpg.max <- input$maxMpgSlider
            #   disp.max <- input$maxDispSlider
            #   tbl.sub <- subset(private$tbl, mpg <= mpg.max & disp <= disp.max)
            #   printf("tbl.sub, %d rows, %d cals", nrow(tbl.sub), ncol(tbl.sub))
            #   private$dtw$setTable(tbl.sub)
            #   })

            #observe({
            #   row.names <- private$dtw$tableSelection()
               #row.names <- rownames(private$tbl)[row.numbers]
               #print(row.names)
            #   private$msgBox$setText(paste(row.names, collapse=", "))
            #   })
            } # server

       ) # public
    ) # class
#--------------------------------------------------------------------------------
# needs this first
deploy <- function()
{
   repos <- options("repos")[[1]]
   stopifnot(sort(names(repos)) == c("BioCann", "BioCsoft", "CRAN"))
   stopifnot(repos$BioCann=="https://bioconductor.org/packages/3.12/data/annotation")
   stopifnot(repos$BioCsoft=="https://bioconductor.org/packages/3.12/bioc")
   stopifnot(repos$CRAN=="https://cran.microsoft.com")

   require(devtools)
   install_github("PriceLab/BioShiny/DataTableWidget", force=TRUE)

   require(rsconnect)

   deployApp(account="hoodlab",
              appName="GATA2-clip-seq",
              appTitle="GATA2-clip-seq",
              appFiles=c("rbpQuickLook.R", "tbl.gata2.RData"),
              appPrimaryDoc="rbpQuickLook.R",
              forceUpdate=TRUE
              )


} # deploy
#--------------------------------------------------------------------------------
app <- DataTableDemoApp$new()
if(grepl("hagfish", Sys.info()[["nodename"]]) & !interactive()){
   runApp(shinyApp(app$ui(), app$server), port=1113)
   } else {
   shinyApp(app$ui(), app$server)
   }


