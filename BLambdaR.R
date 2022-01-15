# Michael Anenburg
# BLambdarR
# See https://lambdar.rses.anu.edu.au for details
# 03 March 2021
# Citable reference: https://doi.org/10.1180/mgm.2020.70

library(shiny)
library(shinyBS)
library(rhandsontable)

# Upload file limit
options(shiny.maxRequestSize = 30*1024^2)

# dput(brewer.pal(5, "Dark2"))
pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")

# default data to show
def.data <- structure(list(La = c(1.447, 12.462, 2.239, 3.479),
                           Ce = c(5.241, 26.642, 5.848, 10.658),
                           Pr = c(1.014, 3.437, 0.947, 1.795),
                           Nd = c(5.722, 14.523, 5.155, 8.99),
                           Sm = c(2.277, 3.464, 1.856, 3.093),
                           Eu = c(1.005, 1.189, 0.742, 1.183),
                           Gd = c(3.285, 3.652, 2.684, 3.612),
                           Tb = c(0.646, 0.598, 0.52, 0.686),
                           Dy = c(4.336, 3.722, 3.416, 4.215),
                           Ho = c(0.906, 0.766, 0.749, 0.92),
                           Er = c(2.715, 2.155, 2.235, 2.596),
                           Tm = c(0.412, 0.325, 0.347, 0.387),
                           Yb = c(2.501, 2.125, 2.282, 2.532),
                           Lu = c(0.351, 0.294, 0.327, 0.358)),
                      class = "data.frame", row.names = c(NA, -4L))

l.header <- c("<em>&lambda;</em><sub>0</sub>",
              "<em>&lambda;</em><sub>1</sub>",
              "<em>&lambda;</em><sub>2</sub>",
              "<em>&lambda;</em><sub>3</sub>",
              "<em>&lambda;</em><sub>4</sub>")
t.header <- c("<em>&tau;</em><sub>1</sub>",
              "<em>&tau;</em><sub>2</sub>",
              "<em>&tau;</em><sub>3</sub>",
              "<em>&tau;</em><sub>4</sub>")

grid.col <- "grey97"
grid.lty <- "solid"

# From O'Neill 2016
radii <- c(1.160, # La
           1.143, # Ce
           1.126, # Pr
           1.109, # Nd
           1.079, # Sm
           1.066, # Eu
           1.053, # Gd
           1.040, # Tb
           1.027, # Dy
           1.015, # Ho
           1.004, # Er
           0.994, # Tm
           0.985, # Yb
           0.977) # Lu

ON16 <- c(0.2472, 0.6308, 0.0950, # La Ce Pr
          0.4793, 0.1542, 0.0592, # Nd Sm Eu
          0.2059, 0.0375, 0.2540, # Gd Tb Dy
          0.0554 ,0.1645, 0.0258, # Ho Er Tm
          0.1684, 0.0251) # Yb Lu

names(radii) <- c("La",
                  "Ce",
                  "Pr",
                  "Nd",
                  "Sm",
                  "Eu",
                  "Gd",
                  "Tb",
                  "Dy",
                  "Ho",
                  "Er",
                  "Tm",
                  "Yb",
                  "Lu")

# Tetrad effect functions, tetrad() is the general structure
tetrad <- function(x) {
  (1-x^2)/2 + sqrt(((1-x^2)^2)/4)
}
# tf() determines size of parabola, 
tf <- function(x, t) { # t is 1, 2, 3, or 4
  x0 <- 55.25 + 3.5*t
  tetrad((x-x0)/1.75)
}
# REE atomic numbers, exclude Pm
z <- (57:71)[-5]
# polynomial fit of radii and atomic numbers
fit <- lm(z ~ poly(radii, 7))
# transformation function that takes ionic radii and outputs atomic number
# note that here it might be simpler to make a simpler function that maps the radius to atomic number,
# e.g. r(La) to 57. However, BLambdaR requires this in proper function form for the least-squares fitting,
# so this form is retained here for consistency.
trn <- function(r) {
  predict(fit, newdata = data.frame(radii = r))
}

bilist <- list("lambda0", "lambda1", "lambda2", "lambda3", "lambda4",
               "tau1", "tau2", "tau3", "tau4",
               "adj. r-squared", "red. chi-sq")

# user interface ----------------------------------------------------------

ui <- fluidPage(
  titlePanel(HTML("BLambdaR: REE patterns to shape coefficients")),
  "Citable reference: ",
  tags$a(href="https://doi.org/10.1180/mgm.2020.70", "Anenburg (2020)."), tags$br(),
  tags$a(href="https://lambdar.rses.anu.edu.au/codeb/instructions.html",
         tags$strong("INSTRUCTIONS - READ THIS FIRST")), tags$br(),
  "Please contact me regarding any issues: ",
  tags$a(href = "mailto:michael.anenburg@anu.edu.au", "michael.anenburg@anu.edu.au"),
  tags$p("Recorded videos of an online workshop on the method are ",
         tags$a(href = "https://www.minersoc.org/ree-software-workshop.html", "available at the Mineralogical Society website", target = "_blank", .noWS = "outside"),
         "."),
  tags$a(href = "http://rses.anu.edu.au", title = "Research School of Earth Sciences", target="_blank",
         tags$div(style = "position:absolute; top:0; right:0",
                  tags$img(src = "https://lambdar.rses.anu.edu.au/codea/ANU_LOGO_CMYK_56mm.jpg",
                           style = "width:215px; height:auto;")
         )
  ),
  fluidRow(
    column(6,
           tabsetPanel(
             tabPanel("Settings",
                      column(6,
                             wellPanel(
                               HTML("A comma separated file (csv) containing 14 columns for each of the REE, in ppm. Each pattern should be in its own row. See <a href=\"https://lambdar.rses.anu.edu.au/codeb/hon.csv\">sample file (the O'Neill 2016 database)</a>."),
                               fileInput("file1", label = "Choose csv file",
                                         multiple = FALSE,
                                         accept = c("text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv")
                               ),
                               "Anomalous elements to exclude from fit:",
                               tags$div(class = "row",
                                        tags$div(class = "col-sm-3",
                                                 bsTooltip(id = "cean", placement = "top",
                                                           title = "When checked, &lsquo;Ce_MFR&rsquo; shows the Ce anomaly: Ce/Ce<sup>*</sup><sub><em>&lambda;&tau;</em></sub>"),
                                                 checkboxInput("cean", label = "Ce", value = FALSE)),
                                        tags$div(class = "col-sm-3",
                                                 bsTooltip(id = "euan", placement = "top",
                                                           title = "When checked, &lsquo;Eu_MFR&rsquo; shows the Eu anomaly: Eu/Eu<sup>*</sup><sub><em>&lambda;&tau;</em></sub>"),
                                                 checkboxInput("euan", label = "Eu", value = TRUE)),
                                        tags$div(class = "col-sm-3",
                                                 bsTooltip(id = "gdan", placement = "top",
                                                           title = "When checked, &lsquo;Gd_MFR&rsquo; shows the Gd anomaly: Gd/Gd<sup>*</sup><sub><em>&lambda;&tau;</em></sub>"),
                                                 checkboxInput("gdan", label = "Gd", value = FALSE)),
                               ),
                               tags$div(class = "row",
                                        tags$div(class = "col-sm-6",
                                                 bsTooltip(id = "norml", placement = "top",
                                                           title = "CI normalisation values from O&rsquo;Neill (2016), Table 1"),
                                                 checkboxInput("norml", label = "Normalise to chondrite", value = TRUE)
                                        ),
                                        tags$div(class = "col-sm-4",
                                                 checkboxInput("fitl4", label = HTML("Fit <em>&lambda;</em><sub>4</sub>"), value = TRUE)
                                        )
                               ),
                               tags$div(class = "row",
                                        tags$div(class = "col-sm-3",
                                                 checkboxInput("fitt1", label = HTML("Fit <em>&tau;</em><sub>1</sub>"), value = FALSE)
                                        ),
                                        tags$div(class = "col-sm-3",
                                                 checkboxInput("fitt2", label = HTML("Fit <em>&tau;</em><sub>2</sub>"), value = FALSE)
                                        ),
                                        tags$div(class = "col-sm-3",
                                                 checkboxInput("fitt3", label = HTML("Fit <em>&tau;</em><sub>3</sub>"), value = FALSE)
                                        ),
                                        tags$div(class = "col-sm-3",
                                                 checkboxInput("fitt4", label = HTML("Fit <em>&tau;</em><sub>4</sub>"), value = FALSE)
                                        )
                               )
                             )
                      ),
                      column(6,
                             wellPanel(
                               tags$strong("Plot:"),
                               tags$div(class = "row",
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plotdata", label = "Data", value = TRUE)),
                                        tags$div(class = "col-sm-4",
                                                 checkboxInput("plotpat", label = "Pattern", value = TRUE)),
                                        tags$div(class = "col-sm-6",
                                                 checkboxInput("plotfit", label = "Fitted curve", value = TRUE)),
                               ),
                               tags$div(class = "row",
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plotl0", label = HTML("<em>&lambda;</em><sub>0</sub>"), value = FALSE)),
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plotl1", label = HTML("<em>&lambda;</em><sub>1</sub>"), value = FALSE)),
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plotl2", label = HTML("<em>&lambda;</em><sub>2</sub>"), value = FALSE)),
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plotl3", label = HTML("<em>&lambda;</em><sub>3</sub>"), value = FALSE)),
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plotl4", label = HTML("<em>&lambda;</em><sub>4</sub>"), value = FALSE)),
                               ),
                               tags$div(class = "row",
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plott1", label = HTML("<em>&tau;</em><sub>1</sub>"), value = FALSE)),
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plott2", label = HTML("<em>&tau;</em><sub>2</sub>"), value = FALSE)),
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plott3", label = HTML("<em>&tau;</em><sub>3</sub>"), value = FALSE)),
                                        tags$div(class = "col-sm-2",
                                                 checkboxInput("plott4", label = HTML("<em>&tau;</em><sub>4</sub>"), value = FALSE)),
                               ),
                               tags$div(class = "row",
                                        tags$div(class = "col-sm-6",
                                                 bsTooltip(id = "params", placement = "top",
                                                           title = "Use &lsquo;legacy&rsquo; only when perfect compatibility with values generated from the spreadsheet in O&rsquo;Neill (2016) are required.<br>For all other cases, use the default &lsquo;full&rsquo;."),
                                                 selectInput("params", "F matrix parameterisation",
                                                             choices = list("full", "legacy"))              
                                        ),
                                        tags$div(class = "col-sm-6",
                                                 bsTooltip(id = "uncert", placement = "top",
                                                           title = "Estimated uncertainty for individual ln(REEs) in percent. See supplementary information in O&rsquo;Neill (2016) for details.<br />Use alongside the &lsquo;&chi;<sup>2</sup><sub><em>&nu;</em></sub> distribution&rsquo; tab."),
                                                 numericInput("uncert", "s(ln[REE])* (%)",
                                                              value = 2, min = 0, step = 0.1))
                               ),
                               tags$div(class = "row",
                                        tags$div(class = "col-sm-6",
                                                 selectInput("bix", "X axis",
                                                             choices = bilist)
                                        ),
                                        tags$div(class = "col-sm-6",
                                                 selectInput("biy", "y axis",
                                                             choices = bilist)
                                        )
                               )
                             )
                      )
             ),
             tabPanel("Change log",
                      fluidRow(wellPanel("2021-03-03: Added link to online workshop.", tags$br(),
                                         "2021-02-19: Added option to plot the measured pattern, on by default.", tags$br(),
                                         "2021-02-18: Plotted fit now takes into account selected anomalies.", tags$br(),
                                         "2021-01-25: Fixed bug that everything was fitted twice, and added a progress indicator.", tags$br(),
                                         "2021-01-24: Updated anomaly notation. Added source code link.", tags$br(),
                                         "2021-01-11: Clarified adjusted r-squared notation.", tags$br(),
                                         "2021-01-09: Added instructions. Updated sample file to use O'Neills database. Fixed bug with missing chi-sq values.", tags$br(),
                                         "2021-01-08: Fixed bug with negative values, treated as missing data.", tags$br(),
                                         "2021-01-03: Fixed bug with zero values.", tags$br(),
                                         "2021-01-02: Tetrads now use corrected parabolas according to polynomial fit to ionic radii", tags$br(),
                                         "2020-12-29: Added option for anomalous Gd.", tags$br(),
                                         "2020-12-23: More decimal digits, improves precision for ppb-level data.", tags$br(),
                                         "2020-12-15: Added bivariate plots.", tags$br(),
                                         "2020-12-14: Added reduced chi-sq tab with histogram, observed, and theoretical distribution.", tags$br(),
                                         "2020-12-11: Now using tau symbols for tetrad coefficients.", tags$br(),
                                         "2020-11-30: Cosmetic upgrades, more detailed file export.", tags$br(),
                                         "2020-11-29: Fixed a small mistake in the input.", tags$br(),
                                         "2020-11-27: Implemented tetrads, added stats: r-squared, se, p-values.", tags$br(),
                                         "2020-11-18: Added option to prevent chondrite normalisation (e.g. when using partition coefficients).", tags$br(),
                                         "2020-11-17: Added MSWD.", tags$br(),
                                         HTML("2020-11-16: Improved data table appearance, added option not to fit <em>&lambda;</em><sub>4</sub>."), tags$br(),
                                         "2020-11-15: Enabled downloads.", tags$br(),
                                         "2020-11-14: Added file upload.", tags$br(),
                                         "2020-11-13: Basic functionality working.", tags$br(),
                                         "2020-11-05: Initial release."))),
             tabPanel(
               "License and source code",
               fluidRow(
                 wellPanel(
                   tags$h2("License"),
                   tags$h4("The MIT License"),
                   tags$p(HTML("Copyright &#169; 2021 Michael Anenburg")),
                   tags$p(HTML("Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the &lsquo;Software&rsquo;), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:")),
                   tags$p("The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software."),
                   tags$p(HTML("THE SOFTWARE IS PROVIDED &lsquo;AS IS&rsquo;, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.")),
                   tags$h2("Source code"),
                   tags$p("Download in R format:"),
                   tags$a("app.R", href="../codeb/app.R"))
               )
             )
           )
    ),
    column(6,
           tabsetPanel(
             tabPanel("REE pattern",
                      fluidRow(
                        column(11,plotOutput(outputId = "reeplot")),
                        column(1,
                               downloadButton(outputId = "patdlpdf", label = "PDF"),
                               downloadButton(outputId = "patdlpng", label = "PNG")
                        )
                      )
             ),
             tabPanel("Bivariate plots",
                      fluidRow(
                        column(11,plotOutput(outputId = "xyplot")),
                        column(1,
                               downloadButton(outputId = "xydlpdf", label = "PDF"),
                               downloadButton(outputId = "xydlpng", label = "PNG")
                        )
                      )
             ),
             tabPanel(HTML("&chi;<sup>2</sup><sub><em>&nu;</em></sub> distribution"),
                      fluidRow(
                        column(11,
                               bsTooltip(id = "chiplot", title = "This is the reduced chi-s distribution (&chi;<sup>2</sup><sub><em>&nu;</em></sub>, or the &chi;<sup>2</sup><sub><em>&nu;</em></sub> divided by the degrees of freedom) of the data. For large datasets, the overall shape of the data distribution (blue) is expected to be roughly similar to the theoretical distribution (red). The position of the theoretical maximum can be adjusted by changing the value of s(ln[REE]), until the distributions match. The resulting value is then an estimate for the total uncertainty in measured REE contents.", placement = "bottom"),
                               plotOutput(outputId = "chiplot")),
                        column(1,
                               downloadButton(outputId = "chidlpdf", label = "PDF"),
                               downloadButton(outputId = "chidlpng", label = "PNG")
                        )
                      )
             )
           )
    )
  ),
  fluidRow(
    wellPanel(
      downloadButton(outputId = "dlcsv", label = "Download results"),
      "MP: measured concentration (ppm), MFR: measured/fit ratio.",
      tags$hr(),
      rHandsontableOutput("hot")
    )
  ),
  tags$div(style="text-align: center",
           tags$a(href="https://biblehub.com/psalms/102-14.htm",
                  tags$img(src="../codea/psalms.svg"),
                  target = "_blank"),
           tags$hr(style = "width:50%"))
)

# server function ---------------------------------------------------------

server <- function(input, output, session) {
  
  curr.row <- reactiveVal()
  
  normval <- reactiveValues()
  
  f.matrix <- reactiveValues()
  
  dfs <- reactiveVal() # degrees of freedom
  
  observe({
    if (input$norml) {
      normval$val <- ON16 # use O'Neill's values for normalisation
    } else {
      normval$val <- rep.int(1, 14) # use "1" for normalisation values
    }
  }
  )
  
  lambda2pat <- function(l0, l1, l2, l3, l4 = 0,
                         t1 = 0, t2 = 0, t3 = 0, t4 = 0) {
    
    # if parameters are not fitted, then the call from reeplot.var() gives NAs.
    # convert them to zeros.
    if (is.na(l4)) {
      l4 <- 0
    }
    if (is.na(t1)) {
      t1 <- 0
    }
    if (is.na(t2)) {
      t2 <- 0
    }
    if (is.na(t3)) {
      t3 <- 0
    }
    if (is.na(t4)) {
      t4 <- 0
    }
    
    # calculate patterns from lambdas
    calc.mat <- apply(t(t(f.matrix[["val"]]) * c(l1,
                                                 l2,
                                                 l3,
                                                 l4)),
                      MARGIN = 1, FUN = sum)
    
    mat.l0 <- calc.mat + l0
    # add tetrads
    print(t1)
    #print(tf(trn(radii), 1))
    tet <- mat.l0 + t1*tf(trn(radii), 1) + t2*tf(trn(radii), 2) + t3*tf(trn(radii), 3) + t4*tf(trn(radii), 4)
    #print(tet)
    return(exp(tet))
  }
  
  
  prmts <- reactive({ # calculate f matrix and orthogonal parameters
    
    if (input$params == "legacy") { # O'Neill's original formulation excluded the radius of Eu
      radii.c <- radii[-6]
    } else {
      radii.c <- radii 
    }
    
    # See appendix in paper for a description of the following section
    N <- length(radii.c)
    a <- sum(radii.c)
    b <- sum(radii.c^2)
    c <- sum(radii.c^3)
    d <- sum(radii.c^4)
    e <- sum(radii.c^5)
    f <- sum(radii.c^6)
    g <- sum(radii.c^7)
    
    # beta
    
    bet <- mean(radii.c)
    
    # gammas
    
    M <- matrix(data = c(-a, N,
                         -b, a),
                nrow = 2, ncol = 2, byrow = TRUE)
    v <- matrix(data = c(-b,
                         -c),
                nrow = 2, ncol = 1, byrow = FALSE)
    vieta <- solve(M, v)
    gam <- Re(polyroot(c(vieta[2],
                         -vieta[1],
                         1)))
    
    # deltas
    
    M <- matrix(data = c(-b, a, -N,
                         -c, b, -a,
                         -d, c, -b),
                nrow = 3, ncol = 3, byrow = TRUE)
    v <- matrix(data = c(-c,
                         -d,
                         -e),
                nrow = 3, ncol = 1, byrow = FALSE)
    vieta <- solve(M, v)
    del <- Re(polyroot(c(-vieta[3],
                         vieta[2],
                         -vieta[1],
                         1)))
    
    # epsilons
    
    M <- matrix(data = c(-c, b, -a, N,
                         -d, c, -b, a,
                         -e, d, -c, b,
                         -f, e, -d, c),
                nrow = 4, ncol = 4, byrow = TRUE)
    v <- matrix(data = c(-d,
                         -e,
                         -f,
                         -g),
                nrow = 4, ncol = 1, byrow = FALSE)
    vieta <- solve(M, v)
    eps <- Re(polyroot(c(vieta[4],
                         -vieta[3],
                         vieta[2],
                         -vieta[1],
                         1)))
    
    # orthogonal terms (f)
    f1 <- radii - bet
    f2 <- (radii - gam[1])*(radii - gam[2])
    f3 <- (radii - del[1])*(radii - del[2])*(radii - del[3])
    f4 <- (radii - eps[1])*(radii - eps[2])*(radii - eps[3])*(radii - eps[4])
    
    # save the f matrix values
    f.matrix[["val"]] <- matrix(data = c(f1, f2, f3, f4),
                                ncol = 4)
    
    # return the orthogonal parameters
    list(bet = bet, gam = gam, del = del, eps = eps)
  })
  
  fit.ree <- function(k) { # least squares fitting done here
    
    ppms <- k # make a copy of the ppm input data
    k <- log(k / normval$val) # normalise
    rad.fit <- radii # make a copy of the ionic radii
    
    # exclude elements from fit, from heavy to light
    if (input$gdan) {
      k <- k[-7]
      rad.fit <- rad.fit[-7]
    }
    if (input$euan) {
      k <- k[-6]
      rad.fit <- rad.fit[-6]
    }
    if (input$cean) {
      k <- k[-2]
      rad.fit <- rad.fit[-2]
    }
    
    
    # construct the orthogonal component functions
    l1f <- function(x) {
      x - prmts()$bet
    }
    l2f <- function(x) {
      (x - prmts()$gam[1])*(x - prmts()$gam[2])
    }
    l3f <- function(x) {
      (x - prmts()$del[1])*(x - prmts()$del[2])*(x - prmts()$del[3])
    }
    
    # construct the function for fitting
    fmla.str <- "k ~ l1f(rad.fit) + l2f(rad.fit) + l3f(rad.fit)"
    
    # add additional functions according to the checked boxes
    if (input$fitl4) {
      l4f <- function(x) {
        (x - prmts()$eps[1])*(x - prmts()$eps[2])*(x - prmts()$eps[3])*(x - prmts()$eps[4])
      }
      fmla.str <- paste0(fmla.str,
                         " + l4f(rad.fit)")
    }
    if (input$fitt1) {
      fmla.str <- paste0(fmla.str,
                         " + tf(trn(rad.fit), 1)")
    }
    if (input$fitt2) {
      fmla.str <- paste0(fmla.str,
                         " + tf(trn(rad.fit), 2)")
    }
    if (input$fitt3) {
      fmla.str <- paste0(fmla.str,
                         " + tf(trn(rad.fit), 3)")
    }
    if (input$fitt4) {
      fmla.str <- paste0(fmla.str,
                         " + tf(trn(rad.fit), 4)")
    }
    fmla <- as.formula(fmla.str)
    
    # least squares fitting!
    fit.mod <- lm(fmla)
    
    # get the various statistics
    summ <- summary(fit.mod)
    
    result <- coef(fit.mod)[1:4] # coefficients
    sserr <- summ$coefficients[, 2] # standard errors
    spvals <- summ$coefficients[, 4] # p values
    
    # make a copy, extract lambda0 to lambda3
    serr <- sserr[1:4]
    pvals <- spvals[1:4]
    
    # extract the others
    if (input$fitl4) {
      result <- c(result, coef(fit.mod)[["l4f(rad.fit)"]])
      serr <- c(serr, sserr[["l4f(rad.fit)"]])
      pvals <- c(pvals, spvals[["l4f(rad.fit)"]])
    } else {
      result <- c(result, NA)
      serr <- c(serr, NA)
      pvals <- c(pvals, NA)
    }
    if (input$fitt1) {
      result <- c(result, coef(fit.mod)[["tf(trn(rad.fit), 1)"]])
      serr <- c(serr, sserr[["tf(trn(rad.fit), 1)"]])
      pvals <- c(pvals, spvals[["tf(trn(rad.fit), 1)"]])
    } else {
      result <- c(result, NA)
      serr <- c(serr, NA)
      pvals <- c(pvals, NA)
    }
    if (input$fitt2) {
      result <- c(result, coef(fit.mod)[["tf(trn(rad.fit), 2)"]])
      serr <- c(serr, sserr[["tf(trn(rad.fit), 2)"]])
      pvals <- c(pvals, spvals[["tf(trn(rad.fit), 2)"]])
    } else {
      result <- c(result, NA)
      serr <- c(serr, NA)
      pvals <- c(pvals, NA)
    }
    if (input$fitt3) {
      result <- c(result, coef(fit.mod)[["tf(trn(rad.fit), 3)"]])
      serr <- c(serr, sserr[["tf(trn(rad.fit), 3)"]])
      pvals <- c(pvals, spvals[["tf(trn(rad.fit), 3)"]])
    } else {
      result <- c(result, NA)
      serr <- c(serr, NA)
      pvals <- c(pvals, NA)
    }
    if (input$fitt4) {
      result <- c(result, coef(fit.mod)[["tf(trn(rad.fit), 4)"]])
      serr <- c(serr, sserr[["tf(trn(rad.fit), 4)"]])
      pvals <- c(pvals, spvals[["tf(trn(rad.fit), 4)"]])
    } else {
      result <- c(result, NA)
      serr <- c(serr, NA)
      pvals <- c(pvals, NA)
    }
    
    # generate the smooth function from fitted coefficients
    fit.pat <- lambda2pat(result[1],
                          result[2],
                          result[3],
                          result[4],
                          result[5],
                          result[6],
                          result[7],
                          result[8],
                          result[9])
    
    fit.ratio <- ppms / (fit.pat * normval$val) # measured to fitted
    names(fit.ratio) <- paste(names(radii), "MFR", sep = "_")
    
    serr <- serr / abs(result) * 100 # calculate se percent
    ssq <- sum(summ$residuals^2) # sum of squares
    red.chi.sq <- ssq / (0.01*input$uncert)^2 / summ$df[2] # reduced chi squared
    
    result <- c(summ$adj.r.squared, red.chi.sq, result, serr, pvals, fit.ratio)
    dfs(summ$df[2]) # save degrees of freedom
    
    isolate({
      curr.row(curr.row() + 1)
      showNotification(ui = paste0("Row ", curr.row(), " complete"),
                       id = "notif", duration = 2, closeButton = FALSE)
    })
    print(result)
    
    return(result)
  }
  
  ree.data <- reactiveValues() # this stores everything, and displayed in the table
  
  # observe -----------------------------------------------------------------
  
  proc.data <- function(temp.ree.data) {
    curr.row(0) # initialise progress counter
    # takes the ppm values and populates all other columns
    temp.ree.data <- cbind(temp.ree.data,
                           t(apply(X = temp.ree.data[, 1:14], 
                                   MARGIN = 1, FUN = fit.ree)))
    colnames(temp.ree.data)[1:14] <- paste(names(radii), "MP", sep = "_")
    colnames(temp.ree.data)[15] <- "adj. <em>r</em><sup>2</sup>"
    colnames(temp.ree.data)[16] <- "&chi;<sup>2</sup><sub><em>&nu;</em></sub>"
    colnames(temp.ree.data)[17:21] <- l.header
    colnames(temp.ree.data)[22:25] <- t.header
    colnames(temp.ree.data)[26:30] <- paste0(l.header, ", se (%)")
    colnames(temp.ree.data)[31:34] <- paste0(t.header, ", se (%)")
    colnames(temp.ree.data)[35:39] <- paste0(l.header, ", p-val")
    colnames(temp.ree.data)[40:43] <- paste0(t.header, ", p-val")
    # save everything
    ree.data[["tableval"]] <- temp.ree.data
  }
  
  observe({ # watch the table
    if (is.null(input$hot)) { # if table does not exist yet
      proc.data(def.data) # use the default data
    }
  })
  
  observeEvent(input$hot$changes$changes, # watch for manual changes in table
               { # take the new data typed in the table 
                 proc.data(hot_to_r(input$hot)[, 1:14])
               })
  
  # watch ticked boxes
  observeEvent(input$fitt1 | input$fitt2 | input$fitt3 | input$fitt4,
               {
                 if (!is.null(input$hot)) {
                   proc.data(hot_to_r(input$hot)[, 1:14])
                 }
               })
  observeEvent(input$fitl4 | input$cean | input$euan | input$gdan | input$norml,
               {
                 if (!is.null(input$hot)) {
                   proc.data(hot_to_r(input$hot)[, 1:14])
                 }
               })
  
  # draw table --------------------------------------------------------------
  
  
  output$hot <- renderRHandsontable({
    tableh <- min(nrow(ree.data[["tableval"]]) * 40 + 20, 400)
    outtable <- rhandsontable(ree.data[["tableval"]], selectCallback = TRUE, height = tableh, digits = 7) # digits important for precision
    outtable <- hot_col(outtable, col = 15:ncol(ree.data[["tableval"]]), readOnly = TRUE)
    outtable <- hot_table(outtable, highlightCol = TRUE, highlightRow = TRUE)
    outtable <- hot_cols(outtable, renderer = "
    function (instance, td, row, col, prop, value, cellProperties) {
    Handsontable.renderers.TextRenderer.apply(this, arguments);
    td.style.color = 'black';

    if (col >= 0 & col <= 13) td.style.background = '#FEFBEB';
    if (col == 14 | col == 15) td.style.background = '#F4EAF7';
    if (col >= 16 & col <= 20) td.style.background = '#FEEFF1';
    if (col >= 21 & col <= 24) td.style.background = '#EEF9F6';
    if (col >= 25 & col <= 29) td.style.background = '#FEEFF1';
    if (col >= 30 & col <= 33) td.style.background = '#EEF9F6';
    if (col >= 34 & col <= 38) td.style.background = '#FEEFF1';
    if (col >= 39 & col <= 42) td.style.background = '#EEF9F6';
    if (col >= 43) td.style.background = '#E7F5FA';

    if (col >= 43 & (value >= 1.05 | value <= 0.95)) td.style.color = 'red';
    if (col >= 43 & (value < 1.05 & value > 1.025)) td.style.color = 'darkred';
    if (col >= 43 & (value > 0.95 & value < 0.975)) td.style.color = 'darkred';
            }")
    #outtable <- hot_col(outtable, col = 1, format="0%")
    outtable <- hot_table(outtable, height = 200)
  })
  
  reeplot.var <- function() {
    if (is.null(input$hot_select$select$r)) {
      plot(1,1,type="n", xaxt = "n", yaxt="n", bty = "n", xlab = "", ylab = "")
      text(1,1, "Select a row", cex = 2)
    } else {
      if (input$hot_select$select$r == input$hot_select$select$r2) {
        
        par(mar = c(2.8, 2.8, 0.1, 0.1),
            mgp = c(1.7, 0.7, 0))
        
        ppms <- ree.data[["tableval"]][input$hot_select$select$r, 1:14]
        lambdas <- unlist(ree.data[["tableval"]][input$hot_select$select$r, 17:21])
        tetrads <- unlist(ree.data[["tableval"]][input$hot_select$select$r, 22:25])
        norm <- ppms / normval$val
        print(tetrads)
        fit.pat <- lambda2pat(l0 = lambdas[1],
                              l1 = lambdas[2],
                              l2 = lambdas[3],
                              l3 = lambdas[4],
                              l4 = lambdas[5],
                              t1 = tetrads[1],
                              t2 = tetrads[2],
                              t3 = tetrads[3],
                              t4 = tetrads[4])
        
        print(fit.pat)
        
        maxs <- max(norm, fit.pat, na.rm = TRUE)*1.1
        mins <- min(norm, fit.pat, na.rm = TRUE)/1.1
        
        plot(1,1, type = "n", bty = "n", log = "y",
             xlim = c(1.16, 0.977), ylim = c(mins, maxs),
             xlab = "Radius / \U00C5", ylab = "CI-normalised (O'N16)")
        
        grid(nx = NA, ny = NULL, lty = grid.lty, col = grid.col) # Horizontal grid
        abline(v = radii, lty = grid.lty, col = grid.col) # Vertical lines at REE
        
        if (input$plotpat) {
          lines(x = radii, y = norm, col = "grey")
        }
        
        if (input$plotdata) {
          points(x = radii, y = norm, cex = 2, pch = 10)
        }
        
        if (input$plotfit) {
          if (input$cean) {
            lines(x = radii[1:3], y = c(fit.pat[1] , norm[2], fit.pat[3]), lwd = 2)
            lines(x = radii[1:3], y = fit.pat[1:3], lwd = 2, lty = "dotted")
          } else {
            lines(x = radii[1:3], y = fit.pat[1:3], lwd = 2)
          }
          
          lines(x = radii[3:5], y = fit.pat[3:5], lwd = 2)
          
          if (input$euan & input$gdan) {
            lines(x = radii[5:8], y = c(fit.pat[5] , norm[6], norm[7], fit.pat[8]), lwd = 2)
            lines(x = radii[5:8], y = fit.pat[5:8], lwd = 2, lty = "dotted")
            lines(x = radii[8:14], y = fit.pat[8:14], lwd = 2)
          } 
          
          if (input$euan & !input$gdan) {
            lines(x = radii[5:7], y = c(fit.pat[5] , norm[6], fit.pat[7]), lwd = 2)
            lines(x = radii[5:7], y = fit.pat[5:7], lwd = 2, lty = "dotted")
            lines(x = radii[7:14], y = fit.pat[7:14], lwd = 2)
          }
          
          if (!input$euan & input$gdan) {
            lines(x = radii[6:8], y = c(fit.pat[6] , norm[7], fit.pat[8]), lwd = 2)
            lines(x = radii[6:8], y = fit.pat[6:8], lwd = 2, lty = "dotted")
            lines(x = radii[5:6], y = fit.pat[5:6], lwd = 2)
            lines(x = radii[8:14], y = fit.pat[8:14], lwd = 2)
          }
          
          if (!input$euan & !input$gdan) {
            lines(x = radii[5:14], y = fit.pat[5:14], lwd = 2)
          }
        }
        
        if (input$plotl0) {
          abline(h = exp(lambdas[1]), col = pal[1], lwd = 1)
        }
        
        if (input$plotl1) {
          lines(x = radii, y = lambda2pat(lambdas[1], lambdas[2], 0, 0, 0, 0, 0, 0, 0), col = pal[2], lwd = 1)
        }
        
        if (input$plotl2) {
          lines(x = radii, y = lambda2pat(lambdas[1], 0, lambdas[3], 0, 0, 0, 0, 0, 0), col = pal[3], lwd = 1)
        }
        
        if (input$plotl3) {
          lines(x = radii, y = lambda2pat(lambdas[1], 0, 0, lambdas[4], 0, 0, 0, 0, 0), col = pal[4], lwd = 1)
        }
        
        if (input$plotl4) {
          if (input$fitl4) {
            l.ty <- "solid"
          } else {
            l.ty <- "dotted"
          }
          lines(x = radii, y = lambda2pat(lambdas[1], 0, 0, 0, lambdas[5], 0, 0, 0, 0), col = pal[5], lwd = 1, lty = l.ty)
        }
        
        if (input$plott1) {
          if (input$fitt1) {
            l.ty <- "solid"
          } else {
            l.ty <- "dotted"
          }
          lines(x = radii[1:4], y = lambda2pat(lambdas[1], 0, 0, 0, 0, tetrads[1], 0, 0, 0)[1:4], col = pal[6], lwd = 1, lty = l.ty)
        }
        if (input$plott2) {
          if (input$fitt2) {
            l.ty <- "solid"
          } else {
            l.ty <- "dotted"
          }
          lines(x = radii[5:7], y = lambda2pat(lambdas[1], 0, 0, 0, 0, 0, tetrads[2], 0, 0)[5:7], col = pal[6], lwd = 1, lty = l.ty)
        }
        if (input$plott3) {
          if (input$fitt3) {
            l.ty <- "solid"
          } else {
            l.ty <- "dotted"
          }
          lines(x = radii[7:10], y = lambda2pat(lambdas[1], 0, 0, 0, 0, 0, 0, tetrads[3], 0)[7:10], col = pal[6], lwd = 1, lty = l.ty)
        }
        if (input$plott4) {
          if (input$fitt4) {
            l.ty <- "solid"
          } else {
            l.ty <- "dotted"
          }
          lines(x = radii[11:14], y = lambda2pat(lambdas[1], 0, 0, 0, 0, 0, 0, 0, tetrads[4])[11:14], col = pal[6], lwd = 1, lty = l.ty)
        }
        
        axis(1, at = radii, labels = names(radii),
             mgp = c(-3, -1.3, 0), tcl = 0.3)
        box()
        
      } else {
        plot(1,1,type="n", xaxt = "n", yaxt="n", bty = "n", xlab = "", ylab = "")
        text(1,1, "Select a single row", cex = 2)
      }
    }
  }
  
  xyplot.var <- function() {
    switch(input$bix,
           "lambda0" = {
             xdata <- ree.data[["tableval"]][, 17]
             xse <- ree.data[["tableval"]][, 17 + 9]
           },
           "lambda1" = {
             xdata <- ree.data[["tableval"]][, 18]
             xse <- ree.data[["tableval"]][, 18 + 9]
           },
           "lambda2" = {
             xdata <- ree.data[["tableval"]][, 19]
             xse <- ree.data[["tableval"]][, 19 + 9]
           },
           "lambda3" = {
             xdata <- ree.data[["tableval"]][, 20]
             xse <- ree.data[["tableval"]][, 20 + 9]
           },
           "lambda4" = {
             xdata <- ree.data[["tableval"]][, 21]
             xse <- ree.data[["tableval"]][, 21 + 9]
           },
           "tau1" = {
             xdata <- ree.data[["tableval"]][, 22]
             xse <- ree.data[["tableval"]][, 22 + 9]
           },
           "tau2" = {
             xdata <- ree.data[["tableval"]][, 23]
             xse <- ree.data[["tableval"]][, 23 + 9]
           },
           "tau3" = {
             xdata <- ree.data[["tableval"]][, 24]
             xse <- ree.data[["tableval"]][, 24 + 9]
           },
           "tau4" = {
             xdata <- ree.data[["tableval"]][, 25]
             xse <- ree.data[["tableval"]][, 25 + 9]
           },
           "adj. r-squared" = {
             xdata <- ree.data[["tableval"]][, 15]
             xse <- NaN
           },
           "red. chi-sq" = {
             xdata <- ree.data[["tableval"]][, 16]
             xse <- NaN
           },
    )
    switch(input$biy,
           "lambda0" = {
             ydata <- ree.data[["tableval"]][, 17]
             yse <- ree.data[["tableval"]][, 17 + 9]
           },
           "lambda1" = {
             ydata <- ree.data[["tableval"]][, 18]
             yse <- ree.data[["tableval"]][, 18 + 9]
           },
           "lambda2" = {
             ydata <- ree.data[["tableval"]][, 19]
             yse <- ree.data[["tableval"]][, 19 + 9]
           },
           "lambda3" = {
             ydata <- ree.data[["tableval"]][, 20]
             yse <- ree.data[["tableval"]][, 20 + 9]
           },
           "lambda4" = {
             ydata <- ree.data[["tableval"]][, 21]
             yse <- ree.data[["tableval"]][, 21 + 9]
           },
           "tau1" = {
             ydata <- ree.data[["tableval"]][, 22]
             yse <- ree.data[["tableval"]][, 22 + 9]
           },
           "tau2" = {
             ydata <- ree.data[["tableval"]][, 23]
             yse <- ree.data[["tableval"]][, 23 + 9]
           },
           "tau3" = {
             ydata <- ree.data[["tableval"]][, 24]
             yse <- ree.data[["tableval"]][, 24 + 9]
           },
           "tau4" = {
             ydata <- ree.data[["tableval"]][, 25]
             yse <- ree.data[["tableval"]][, 25 + 9]
           },
           "adj. r-squared" = {
             ydata <- ree.data[["tableval"]][, 15]
             yse <- NaN
           },
           "red. chi-sq" = {
             ydata <- ree.data[["tableval"]][, 16]
             yse <- NaN
           },
    )
    if (is.na(xdata) | is.na(ydata)) {
      plot(1,1,type="n", xaxt = "n", yaxt="n", bty = "n", xlab = "", ylab = "")
      text(1,1, "Error: Attempting to plot a non-fitted parameter", cex = 2)
    } else {
      plot(x = xdata, y = ydata, type = "n", bty = "n",
           xlab = input$bix, ylab = input$biy)
      grid()
      if (!is.nan(xse)) {
        segments(x0 = xdata - xdata*xse/100, y0 = ydata,
                 x1 = xdata + xdata*xse/100, y1 = ydata,
                 col = "dimgrey")
        segments(y0 = ydata - ydata*yse/100, x0 = xdata,
                 y1 = ydata + ydata*yse/100, x1 = xdata,
                 col = "dimgrey")
      }
      points(x = xdata, y = ydata, pch = 20)
      box()
    }
  }
  
  chiplot.var <- function() {
    brk <- 50
    # na.omit to remove cases where the fitted REEs are so low that r and chi squared are not reported
    chisqs <- na.omit(ree.data[["tableval"]][, 16])
    dens <- density(chisqs)
    histval <- hist(chisqs, breaks = brk)
    xs <- seq(0, max(chisqs)*dfs(), length = 100*dfs())
    chidist <- dchisq(xs, df = dfs())
    maxval <- max(histval$counts)
    hist(chisqs, breaks = brk, main = "", xlab = "Reduced chi-squared")
    densfac <- max(dens$y) / max(maxval)
    distfac <- max(chidist) / max(maxval)
    lines(y = dens$y / densfac, x = dens$x, col = "blue")
    lines(y = chidist / distfac, x = xs / dfs(), col = "red")
    legend(x = "topright", legend = c("observed", "theoretical"),
           title = "distribution", col = c("blue", "red"), lty = 1)
  }
  
  output$reeplot <- renderPlot({
    print(reeplot.var())
  })
  
  output$xyplot <- renderPlot({
    print(xyplot.var())
  })
  
  output$chiplot <- renderPlot({
    print(chiplot.var())
  })
  
  observeEvent(input$file1,{
    uploaded()
  })
  
  # file upload -------------------------------------------------------------
  
  
  uploaded <- reactive({
    
    req(input$file1)
    
    # Expects a csv file in a certain format
    df <- read.csv(input$file1$datapath,
                   header = FALSE,
                   col.names = names(radii),
                   colClasses = rep("numeric", 14)) # only take numbers
    # Treat zeros as missing values
    df[df <= 0] <- NA
    # processes uploaded data
    proc.data(df)
    return(df)
  })
  
  
  # download handlers -------------------------------------------------------
  
  output$patdlpdf <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%M%e%H%M%S"), ".pdf")
    },
    content = function(file) {
      pdf(file)
      reeplot.var()
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  output$patdlpng <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%M%e%H%M%S"), ".png")
    },
    content = function(file) {
      png(file)
      reeplot.var()
      dev.off()
    },
    contentType = "image/png"
  )
  
  output$xydlpdf <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%M%e%H%M%S"), ".pdf")
    },
    content = function(file) {
      pdf(file)
      xyplot.var()
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  output$xydlpng <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%M%e%H%M%S"), ".png")
    },
    content = function(file) {
      png(file)
      xyplot.var()
      dev.off()
    },
    contentType = "image/png"
  )
  
  output$chidlpdf <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%M%e%H%M%S"), ".pdf")
    },
    content = function(file) {
      pdf(file)
      chisq.var()
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  output$chidlpng <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%M%e%H%M%S"), ".png")
    },
    content = function(file) {
      png(file)
      chisq.var()
      dev.off()
    },
    contentType = "image/png"
  )
  
  output$dlcsv <- downloadHandler(
    filename = function() {
      paste0(format(Sys.time(), "%Y%M%e%H%M%S"), ".csv")
    },
    content = function(file) {
      
      cat("BLambdaR output version 2021-03-03\nhttps://lambdar.rses.anu.edu.au\n", file = file)
      
      if (input$cean) {
        cat("Ce considered anomalous: not used in fit\n", file = file, append = TRUE)
      } else {
        cat("Ce not considered anomalous: used in fit\n", file = file, append = TRUE)
      }
      if (input$euan) {
        cat("Eu considered anomalous: not used in fit\n", file = file, append = TRUE)
      } else {
        cat("Eu not considered anomalous:  used in fit\n", file = file, append = TRUE)
      }
      if (input$gdan) {
        cat("Gd considered anomalous: not used in fit\n", file = file, append = TRUE)
      } else {
        cat("Gd not considered anomalous:  used in fit\n", file = file, append = TRUE)
      }
      if (input$params == "legacy") {
        cat("Legacy F-matrix parameterisation\n", file = file, append = TRUE)
      } else {
        cat("Full F-matrix parameterisation\n", file = file, append = TRUE)
      }
      if (input$fitl4) {
        cat("Lambda 4 included in fitting\n", file = file, append = TRUE)
      } else {
        cat("Lambda 4 not included in fitting\n", file = file, append = TRUE)
      }
      if (input$fitt1) {
        cat("Tetrad 1 included in fitting\n", file = file, append = TRUE)
      } else {
        cat("Tetrad 1 not included in fitting\n", file = file, append = TRUE)
      }
      if (input$fitt2) {
        cat("Tetrad 2 included in fitting\n", file = file, append = TRUE)
      } else {
        cat("Tetrad 2 not included in fitting\n", file = file, append = TRUE)
      }
      if (input$fitt3) {
        cat("Tetrad 3 included in fitting\n", file = file, append = TRUE)
      } else {
        cat("Tetrad 3 not included in fitting\n", file = file, append = TRUE)
      }
      if (input$fitt4) {
        cat("Tetrad 4 included in fitting\n", file = file, append = TRUE)
      } else {
        cat("Tetrad 4 not included in fitting\n", file = file, append = TRUE)
      }
      cat(paste0("s(ln[REE])=", input$uncert, "%"), file = file, append = TRUE)
      
      cat("\n\n", file = file, append = TRUE)
      
      dat <- hot_to_r(input$hot)
      colnames(dat)[15] <- "adj_r-sqrd"
      colnames(dat)[16] <- "red_chi-sqrd"
      colnames(dat)[17:21] <- c("lambda0", "lambda1", "lambda2", "lambda3", "lambda4")
      colnames(dat)[22:25] <- c("tau1", "tau2", "tau3", "tau4")
      colnames(dat)[26:30] <- paste0(c("lambda0", "lambda1", "lambda2", "lambda3", "lambda4"), "_se%")
      colnames(dat)[31:34] <- paste0(c("tau1", "tau2", "tau3", "tau4"), "_se%")
      colnames(dat)[35:39] <- paste0(c("lambda0", "lambda1", "lambda2", "lambda3", "lambda4"), "_pval")
      colnames(dat)[40:43] <- paste0(c("tau1", "tau2", "tau3", "tau4"), "_pval")
      
      write.table(dat, file, append = TRUE, sep = ",", row.names = FALSE, na = "")
    }
  )
}

# shinyApp() --------------------------------------------------------------

shinyApp(ui = ui, server = server)