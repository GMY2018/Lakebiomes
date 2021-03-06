
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
library(gridExtra)
library(shiny)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(maps)


mc=c("#0000CD", "#1E90FF","#ADD8E6","#EEC900" ,"#FF8C00","#EE0000", "#EEA9B8","#9ACD32", "#228B22")[order(c(4, 5, 9, 1, 2, 3, 6, 8, 7))]
load("globalpos.Rdata")
load("globalclass.Rdata")
load("efd.Rdata")
load("foy.Rdata")
load("meanefd.Rdata")

load("globalflakemat.Rdata")
load("uncert.Rdata")
load("uncertall.Rdata")
skip <- which(rowSums(globalflakemat==0)>410)
uall <- round(uncertall[-skip, ],2)
load("arcfd.Rdata") # the seasonal fda object for the arclake data
load("arc_fpca.Rdata")
load("myqda.Rdata")
library(MASS)
library(fda)
pos <- data.frame(globalpos)
names(pos) <- c("x","y")
gp <- data.frame(cbind(x=globalpos[,1],y=globalpos[,2], group=globalclass ))[!is.na(globalclass),]
load("omeanfd.Rdata")
gclass <- globalclass[-skip]
post <- rep(NA, length(gclass))
for(i in 1:length(post)){ post[i] <- uall[i, as.numeric(gclass[i])]}

cnames <- c("Northern cold",
            "Northern cool",
            "Northern temperate",
            "Northern warm",
            "Northern hot",
            "Equatorial hot",
            "Southern hot",
            "Southern warm",
            "Southern temperate")

codes <- c("NC", "NO", "NT", "NW", "NH", "EH", "SH", "SW", "ST" )
cols <- c("Dark blue", "Mid blue", "Light blue", "Yellow", "Orange", "Red", "Pink", "Light Green","Dark Green")
dat <- data.frame(cbind("Name"=cnames, "Code"=codes, "Colour"=cols))


sumfun <- function(x){data.frame(round(cbind(min=min(x, na.rm=TRUE), 
                                             max=max(x, na.rm=TRUE), 
                                             mean=mean(x, na.rm=TRUE),
                                             percent.zeros=sum(x==0)*100/length(x)),2))}

lakedat <- do.call(rbind, apply(globalflakemat,1,sumfun))

makeamap <- function(){
  mp <- NULL
  mapWorld <- borders("world", colour="gray50", fill="gray50")
  mp <- ggplot() +   mapWorld + labs(x="Longitude", y="Latitude") + theme(legend.title=element_text(size=rel(1.75)) , legend.text=element_text(size=rel(1.75)),axis.text = element_text(size = rel(1.75)), axis.title = element_text(size = rel(1.75)))
  gp_ <- gp
  gp_$group <- recode(gp_$group, "4"="NC", "5"="NO", "9"="NT",
                      "1"="NW", "2"="NH", "3"="EH",
                      "6"="SH", "8"="SW", "7"="ST") 
  
  mp1 <- mp+ 
    geom_point(data=gp_, aes(x=gp_$x, y=gp_$y, col=factor(gp_$group, levels=c("NW", "NH", "EH", "NC", "NO", "SH", "ST", "SW", "NT"))),
               size=2, shape=15) + 
    labs(color = "Group")  + 
    scale_color_manual(values=mc) + guides(color = guide_legend(override.aes = list(size=5)))
  mp1
  
}


makeuncertmap <- function(){
  mp <- NULL
  mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
  mp <- ggplot() +   mapWorld + labs(x="Longitude", y="Latitude")
  gp2 <- data.frame(cbind(x=globalpos[-skip,1],y=globalpos[-skip,2], post=post ))
  mp2 <- mp+ geom_point(data=gp2, aes(x=gp2$x, y=gp2$y, color=gp2$post) ,size=2, shape=15)+labs(color = "Post")  + scale_colour_gradientn(colours = rev(terrain.colors(10)))
  mp2 + theme(legend.title=element_text(size=rel(1.75)) , legend.text=element_text(size=rel(1.75)),axis.text = element_text(size = rel(1.75)), axis.title = element_text(size = rel(1.75)))
}

rc <- c("x4"="NC","X5"="NO", "X9"="NT", "X1"="NW","X2"="NH",
        "X3"="EH","X6"="SH", "X8"="SW",	"X7"="ST")

library(dplyr)


gcr <- recode(globalclass, "4"="NC","5"="NO", "9"="NT", "1"="NW","2"="NH",
              "3"="EH","6"="SH", "8"="SW",	"7"="ST") 

plotthemean <- function(){
  meanline <- data.frame(cbind(x=as.numeric(rownames(meanefd)),
                               y=meanefd))
  meltdf <- melt(meanline,id="x")
  
  meltdf$variable <- recode(meltdf$variable, "X4"="NC","X5"="NO", "X9"="NT", "X1"="NW","X2"="NH",
                            "X3"="EH","X6"="SH", "X8"="SW",	"X7"="ST") 
  p <- ggplot(meltdf, aes(x=x, y=value)) +labs(color = "Group")  + scale_color_manual(values=mc)
  d <- p + geom_line(data=meltdf,
                     aes(x=x, y=value,colour=variable,
                         group=variable), size=rel(1.5), 
                     inherit.aes = FALSE) +
    labs(x="Week of Year", y="LSWT Celsius")
  d+ theme(legend.title=element_text(size=rel(1.5)) , 
           legend.text=element_text(size=rel(1.5)),
           axis.text = element_text(size = rel(1.1)), 
           axis.title = element_text(size = rel(1.1)))
  
  
  
}


plotdat <- function(i, am=FALSE, addmeans=NULL){
  tc <- data.frame(cbind(y=globalflakemat[i,], x=foy))
  
  lakeinfo <- paste(paste(c("min=", "max=", "mean=", "% zeros="), lakedat[i,]), collapse = " , " )
  ff <- gcr[i]
  fc <- globalclass[i]
  po <- globalpos[i,]
  u <- round(uncertall[i,],2)
  names(u) <- recode(names(u), "4"="NC","5"="NO", "9"="NT", "1"="NW","2"="NH",
                     "3"="EH","6"="SH", "8"="SW",	"7"="ST") 
  nu <- paste(names(u), u, sep=":", collapse = " , ")
  
  # nuu <- paste0(strsplit(nu, ", NO")[[1]][1], "\n NO",strsplit(nu, ", NO")[[1]][2]) 
  datline <- data.frame(cbind(nx=as.numeric(rownames(efd)), ny=efd[,i]))
  p <- ggplot(tc, aes(x, y))+ theme(legend.position="none", legend.title=element_text(size=rel(1.25)) ,
                                    plot.title=element_text(size=rel(1.5), hjust=0.5),
                                    legend.text=element_text(size=rel(1.25)),
                                    axis.text = element_text(size = rel(1.25)),
                                    axis.title = element_text(size = rel(1.25)),
                                    plot.caption=element_text(size=rel(1.4), hjust=0.5, color="black"),
                                    plot.subtitle=element_text(size=rel(1.4), hjust=0.5, color="black")) +
    labs(title=paste("Lon =", po[1],", Lat =",po[2], ", Pred Group=",ff), subtitle=paste("Post:", nu), caption=paste(lakeinfo))
  d <- p + geom_point() + labs(x="Week of Year", y="LSWT Celsius")
  z <- d + geom_line(data=datline, aes(x=datline$nx, y=datline$ny),size=rel(2), col="purple", inherit.aes = FALSE)+ ylim(0, 35)
  
  if (am) {
    meanline <- data.frame(cbind(x=as.numeric(rownames(meanefd)), y=meanefd[,fc]))
    z <- z + geom_line(data=meanline, aes(x=x,y=y), color=mc[fc], size=rel(2), inherit.aes = FALSE)
  }
  z
  
  if (!is.null(addmeans)) {
    iadm <- unlist(paste(addmeans))
    iadm <- as.numeric(recode(iadm, "NC"="4","NO"="5", "NT"="9", "NW"="1","NH"="2",
                              "EH"="3","SH"="6", "SW"="8",	"ST"="7") )
    for(ii in iadm){
      meanline <- data.frame(cbind(x=as.numeric(rownames(meanefd)), y=meanefd[,ii]))
      z <- z + geom_line(data=meanline, aes(x=x,y=y), color=mc[ii], size=rel(2), inherit.aes = FALSE)
    }
  }
  z
  
}





seasonalfd <- function(dat, wk){
  fulllam <- df2lambda(c(wk), basis=create.bspline.basis(c(0,52), nbasis=53), df=12)
  
  rounded <- round(unlist(lapply(split(dat, wk), mean)),2)
  zeros= as.numeric(names(rounded[rounded<0.5]))+1
  shortBasis <- create.bspline.basis(c(0,52), nbasis=53)
  if (length(zeros)==0){
    ff <- Data2fd(y=dat, argvals=wk, basis=shortBasis, lambda=fulllam)
    nco <- ff$coefs}
  
  if (length(zeros)>0){
    mylam <- df2lambda(c(wk), basis=create.bspline.basis(c(0,52), nbasis=53-length(zeros)), df=12)
    redshortBasis <- create.bspline.basis(c(0,52), nbasis=53, dropind=zeros)
    ff <- Data2fd(y=dat, argvals=wk, basis=redshortBasis, lambda=mylam)
    nco <- rep(0, 53)
    nco[-zeros] <- ff$coefs
  }
  
  nco
}
#td<- read.csv("Z://test_one.csv", stringsAsFactors = FALSE)
#indataset=td

gettestdat <- function(indataset,i){
  lname <- colnames(indataset)[i+1]
  date <- strptime(indataset[,1], format="%d/%m/%Y")
  if(any(is.na(date))) {stop("There are missing dates!")}
  wk <- date$yday%/%7 
  # foy <- indataset[,2]
  rem <- is.na(indataset[,i+1])
  newlswt  <- indataset[!rem,i+1]
  newwk    <- wk[!rem]
  testdata <- data.frame(cbind(lswt=newlswt, wk=newwk))
  list(testdata=testdata, lname=lname)
}


predplotting <- function(testdat2, am=FALSE, addmeans=NULL){
  testdat <- testdat2$testdata
  lname <- testdat2$lname
  nco <- seasonalfd(dat=testdat$lswt, wk=testdat$wk)
  fdo <- fd(nco , create.bspline.basis(c(0,52), nbasis=53))
  ## get the projected scores and set as a data frame (necessary for the predict function)
  sc <- t(inprod(arc.fpc$harmonics, fdo-mean(arcfd)))
  colnames(sc) <- c("sc1", "sc2")
  sc <- as.data.frame(sc)
  pred <- predict(myqda, newdat=sc); cc <- pred$class; pc <- round(pred$posterior,2)
  ccn <- paste(recode(cc, "4"="NC","5"="NO", "9"="NT", "1"="NW","2"="NH","3"="EH","6"="SH", "8"="SW",	"7"="ST") )
  
  #plot( testdat$wk, testdat$lswt, ylim=c(0,40), xlab="Week of Year", ylab="LSWT", cex.axis=1.5, cex.lab=1.5)
  tc <- data.frame(cbind(x=testdat$wk,y=testdat$lswt))
  colnames(pc) <-  recode(colnames(pc), "4"="NC","5"="NO", "9"="NT", "1"="NW","2"="NH","3"="EH","6"="SH", "8"="SW",	"7"="ST") 
  tt <- paste0(colnames(pc), ": ", round(pred$post,2), collapse=", ")
  p <- ggplot(tc, aes(x, y))+ ylim(0, 35) +
    theme(legend.position="none", 
          plot.title=element_text(size=rel(1.5), hjust=0.5), 
          axis.text = element_text(size = rel(1)),
          axis.title = element_text(size = rel(1)))  + labs(title=paste(lname, ": predicted class :", ccn, "\n", tt))
  
  d <- p + labs(x="Week of Year", y="LSWT Celsius")
  datline <- data.frame(cbind(nx=seq(0,52,0.5), ny=c(eval.fd(fdo,seq(0,52,0.5)))))
  z <- d + geom_line(data=datline, aes(x=datline$nx, y=datline$ny),size=rel(1.5), col="purple", inherit.aes = FALSE)
  z <- z+ geom_point()
  z
  
  if (am) {
    meanline <- data.frame(cbind(x=as.numeric(rownames(meanefd)), y=meanefd[,as.numeric(cc)]))
    z <- z + geom_line(data=meanline, aes(x=x,y=y), color=mc[as.numeric(cc)], size=rel(1.5), inherit.aes = FALSE)
    z
  }
  
  if (!is.null(addmeans)) {
    iadm <- unlist(paste(addmeans))
    iadm <- as.numeric(recode(iadm, "NC"="4","NO"="5", "NT"="9", "NW"="1","NH"="2",
                              "EH"="3","SH"="6", "SW"="8",	"ST"="7") )
    for(ii in iadm){
      meanline <- data.frame(cbind(x=as.numeric(rownames(meanefd)), y=meanefd[,ii]))
      z <- z + geom_line(data=meanline, aes(x=x,y=y), color=mc[ii], size=rel(2), inherit.aes = FALSE)
    }
  }
  z
}


basicmap<- function(ppp){
  plot(globalpos, col=c(mc)[globalclass], pch=15, xlab="",ylab="", bty="n")
  points(globalpos[as.numeric(ppp), 1],globalpos[as.numeric(ppp),2 ], lwd=1)
  abline(h=globalpos[as.numeric(ppp), 2], lty=2)
  abline(v=globalpos[as.numeric(ppp), 1], lty=2)
}



alltheplots <- function(indataset,am=FALSE, addmeans=NULL){
  nc <- ncol(indataset)-1
  mylist <- NULL
  for(i in 1:nc){
    testdat2 <- gettestdat(indataset, i)
    mylist[[i]] <- predplotting(testdat2,am, addmeans)
  }
  grid.arrange(grobs=mylist,nrow=length(mylist))
}



ui <- fluidPage(
  
  # Application title
  titlePanel("Towards Lake Biomes?"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      
      
      
      
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
      
      
      checkboxInput("checkbox", label = "Show predicted mean", value = FALSE),
      
      selectizeInput('addmean', 
                     label=h6("Add Group Mean Curves"),
                     choices=c("NC","NO", "NT", "NW","NH","EH","SH", "SW",	"ST"), multiple = TRUE,
                     options = list(
                       placeholder = 'Please select',
                       onInitialize = I('function() { this.setValue(""); }'))),
      
      
      
      #dateRangeInput("inDateRange", "Input date range"),
      
      plotOutput('meanplot'),
      tableOutput('to'),
      
      withTags({div(class="header", checked=NA, 
                    h4("This app shows the predicted lake class for 
                       each grid cell (~2 degree) agross the globe. The time 
                       series which were used to obtain these predictions
                       were based on the Flake data produced by Iestyn.
                       The prediction is based on a smooth curve.  
                       Clicking on a point on the map will produce
                       a plot of the Flake data (points) along with the
                       smoothed curve (shown by a solid purple line).
                       There is also the option to add on
                       the predicted group mean for the location selected and any additional group means that might be of interest.
                       
                       'Post' indicates the posterior probability of cluster membership.   
                       
                       Locations where there is a 'flat' data series at zero have not
                       been included in the predictions due to excessive periods of ice cover, 
                       .i.e. there is no predicted class for these locations. If a location on the map 
                       is selected where no data are available, the error 'no applicable method for 
                       recode applied to an object of class NULL' will be displayed"),
                    
                    
                    h4("For lake prediction the data must be in the form of a csv file 
                       with the first column containing sample date in dd/mm/yyyy format 
                       and the subsqeuent columns containing lswt measurements in celsius. 
                       Any missing LSWT values should be denoted NA. An error will be given if any of the dates are missing "),   
                    
                    h5(b("NOTE:"), "The app is a 'beta version' :) "))
        
      })
      ),
    
    
    
    
    mainPanel( 
      tabsetPanel(
        tabPanel("Main",
                 plotOutput("mapplot", click = "plot_click"),
                 plotOutput("plot2")
        ),
        
        tabPanel("Prediction",
                 
                 plotOutput("predplot")
        ),
        
        tabPanel('User Input Location',   
                 
                 fluidRow(
                   column(2, numericInput("lat", "Longitude:",  NULL, min = -180, max = 180)),
                   column(2,numericInput("lon", "Latitude:", NULL, min = -55, max = 82))
                 ),
                 
                 
                 plotOutput("extramap") ,
                 plotOutput("userplot") 
        ), 
        
        
        tabPanel('Maps',   
                 #  plotOutput("mapplot2"),
                 plotOutput("uncertplot"),
                 imageOutput("image")
        )
        
      ))
      ))





server <- function(input, output, session) {
  
  ## get data
  mydat <- reactive({
    inFile <- input$file1
    if (is.null(inFile))  return(NULL)
    read.csv(inFile$datapath, header=TRUE, stringsAsFactors = FALSE)#input$header)
  })
  
  observe({
    date <- strptime(mydat()[,1], format="%d/%m/%Y")
    
    updateDateRangeInput(session, "inDateRange",
                         label = paste("Date range"),
                         start = sort(date)[1],
                         end = sort(date)[length(date)] )
  })
  
  observe({
    validate(need(input$file1, 'please upload data file, missing values should be included as NA'))
    
    output$predplot <- renderPlot({
      validate(need(input$file1, 'please upload data file, missing values should be included as NA'))
      #testdata <- data.frame(gettestdat(mydat()))
      addmeans <- input$addmean
      amm <- input$checkbox
      pr <- alltheplots(data.frame(mydat()),am=amm, addmeans)
    },
    height = 400*(ncol(mydat())-1))
  })
  
  output$to <- renderTable({dat})
  
  output$image <- renderImage({
    width  <- session$clientData$output_image_width
    height <- session$clientData$output_image_height
    filename <- normalizePath("gms.png")
    list(src = filename, width = width,
         height = height*2 ) }, deleteFile = FALSE)
  
  
  output$extramap <- renderPlot({
    validate(need(input$lat, 'please enter a latitude & longitude'))
    validate(need(input$lon, 'please enter a latitude & longitude'))
    pp <- c(input$lat, input$lon)
    ppp <- which.min(sqrt((globalpos[,1]-pp[1])^2+(globalpos[,2]-pp[2])^2))[1]
    basicmap(ppp)
  })
  
  output$mapplot2 <- renderPlot(makeamap())
  
  output$mapplot <- renderPlot(makeamap())
  output$meanplot <- renderPlot(plotthemean())  
  output$uncertplot <- renderPlot(makeuncertmap())
  
  
  output$userplot <- renderPlot({
    validate(need(input$lat, ''))
    validate(need(input$lon, ''))
    
    pp <- c(input$lat, input$lon)
    ppp <- which.min(sqrt((globalpos[,1]-pp[1])^2+(globalpos[,2]-pp[2])^2))[1]
    plotdat(ppp, input$checkbox, input$addmean)
  })
  
  
  
  output$plot2 <- renderPlot({
    validate(need(input$plot_click, 'Select a point from the map'))
    pp <- nearPoints(pos, input$plot_click, maxpoints = 1,xvar = "x", yvar = "y", allRows = TRUE)
    
    ppp <- pp$selected_
    if(is.na(ppp)) cat("There are no nearby points")
    plotdat(ppp, input$checkbox, input$addmean)
    
  })
  
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# Run the application 
shinyApp(ui = ui, server = server)

