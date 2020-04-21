library(shiny)
library(shinydashboard)
library("htmltools")
library("bsplus")

ui <- dashboardPage(
  dashboardHeader(title = 'Coming out of COVID-19', titleWidth = 250),
  dashboardSidebar(disable = 'TRUE'),
  dashboardBody(
    fluidRow(
      box(width = 4,
          title = "Controls",
          bs_accordion(id = "sliderou") %>% 
            bs_append(title =  sliderInput("slidero", "Transmissibility level for within the over-65 population:",
                                           min = 0, max = 3, animate = TRUE, value = 2.8, step = 0.1),
                      content = "Transmission rate within over-65 population: This determine the rate at which encounters between infected
and susceptible individuals lead to new infections. Due to Assumption 3 (Society will eventually return to near-normal transmission levels),
we assume that after 1 year, these transmission rates will return to levels that imply an R0-value of 2.8 in a completely susceptible population.
The R0 value correspond to the median estimate across a several studies.")%>% 
            bs_append(title =  sliderInput("slidery",
                                           "Transmissibility level for within the under-65 population:",
                                           min = 0, max = 3, animate = TRUE, value = 2.8, step = 0.1),
                      content = "Transmission rate within under-65 population: This determine the rate at which encounters between infected
and susceptible individuals lead to new infections. Due to Assumption 3 (Society will eventually return to near-normal transmission levels),
we assume that after 1 year, these transmission rates will return to levels that imply an R0-value of 2.8 in a completely susceptible population.
The R0 value correspond to the median estimate across a several studies.") %>%
            bs_append(title = sliderInput("age.mortality.ratio", 
                                          "Age-relative mortality rate:",
                                          min = 0, max = 0.25, animate = TRUE, value = (1/50), step = 0.02),
                      content = "Assume the mortality rate for our older population is 50 times greater than for our younger population. Thus the rate is 1/50.") %>%
            bs_append(title = sliderInput("N", "Total population size:",
                                          min = 0, max = 300e6, animate = TRUE, value = 300e6, step = 1e6),
                      content = "Total population size: Default is 300 million, roughly the population of the US.") %>%
            bs_append(title = sliderInput("sf", "Susceptible fraction:",
                                          min = 0.01, max = 1, animate = TRUE, value = 0.25, step = 0.01), 
                      content = "Susseptible fraction times total population size is the total number of susceptible individuals at time 0. Susceptible fraction is set to 0.25 based
on observations of the total infected population in flu pandemics. Of
course the actual value may be much closer to 1 however since nonsusceptible individuals do not interact with the other compartments in
the SIR model this has no effect on the dynamics of the dynamics, only
the absolute numbers.") %>%
            bs_append(title = sliderInput("overcrowdthresh", "Hospital overcrowding threshold:",
                                          min = 0, max = 1e6, animate = TRUE, value = 5e5 , step = 1e5),
                      content = "Medical system capacity: The US is estimated to have a total capacity
of 728,000 hospital beds. Of course not all of these beds are available
for COVID-19 patients. To model the effects of overburdened medical
systems, assumption is that above a threshold of 500,000 infected cases, mortality increases by a factor of 2. It is shown in the sensitivity analysis that
the model is not sensitive to these choices.")
      ),
      box(width = 8, title = 'Results',
          plotOutput("plot1", height = 400)),
      HTML("<br><br>"),
      box(width = 8,
          'This  tool is built for interactive visualization of the results of the paper', a("Fighting COVID-19: the heterogeneous transmission thesis.",  target="_blank", href =  "https://math.cmu.edu/~wes/covid.html#footM")
      )
    )
  )
)

server <- function(input, output) {
  source('simFuncs.R')
  output$plot1 <- renderPlot({
    o <- input$slidero
    y <- input$slidery
    oct <- input$overcrowdthresh
    base.mortality.ratio <- input$base.mortality.ratio
    age.mortality.ratio <- input$age.mortality.ratio
    N <- input$N
    sf <- input$sf
    overcrowdthresh <- input$overcrowdthresh
    m <- simulateSIR2pop(bo.mult.init = o, by.mult.init = y,
                         age.mortality.ratio = age.mortality.ratio, 
                         N = N, 
                         sf = sf,
                         overcrowdthresh = overcrowdthresh)
    plotResults(m, step.size = 1, oct = oct)
  })
}

shinyApp(ui, server)