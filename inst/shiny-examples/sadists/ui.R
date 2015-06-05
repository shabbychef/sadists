# Created: 2015.05.18
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

library(shiny)

ui <- shinyUI(fluidPage(
  titlePanel("Simulations"),
  sidebarLayout(
    sidebarPanel(
      h3("parameters"),
      selectInput("distro", "Distribution:", 
									choices=c("dnbeta","dneta","dnf","dnt","kprime","lambap",
														"prodchisqpow","proddnf","sumchisqpow","sumlogchisq","upsilon"),
									selected="upsilon",
									multiple=FALSE),
			conditionalPanel(
				condition = "input.distro == 'dnbeta'",
				numericInput("dnbeta_df1", "df1:", min=1, max=Inf, value=50, step=0.1),
				numericInput("dnbeta_df2", "df2:", min=1, max=Inf, value=100, step=0.1),
				numericInput("dnbeta_ncp1", "ncp1:", min=0, max=Inf, value=1, step=0.001),
				numericInput("dnbeta_ncp2", "ncp2:", min=0, max=Inf, value=2, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'dneta'",
				numericInput("dneta_df", "df:", min=-Inf, max=Inf, value=50, step=0.1),
				numericInput("dneta_ncp1", "ncp1:", min=0, max=Inf, value=1, step=0.001),
				numericInput("dneta_ncp2", "ncp2:", min=0, max=Inf, value=2, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'dnf'",
				numericInput("dnf_df1", "df1:", min=1, max=Inf, value=50, step=0.1),
				numericInput("dnf_df2", "df2:", min=1, max=Inf, value=100, step=0.1),
				numericInput("dnf_ncp1", "ncp1:", min=0, max=Inf, value=1, step=0.001),
				numericInput("dnf_ncp2", "ncp2:", min=0, max=Inf, value=2, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'dnt'",
				numericInput("dnt_df", "df:", min=1, max=Inf, value=50, step=0.1),
				numericInput("dnt_ncp1", "ncp1:", min=-Inf, max=Inf, value=1, step=0.001),
				numericInput("dnt_ncp2", "ncp2:", min=0, max=Inf, value=2, step=0.001)
			),
			conditionalPanel(
				condition = "input.distro == 'kprime'",
				numericInput("kprime_a", "a:", min=0, max=Inf, value=5, step=0.1),
				numericInput("kprime_b", "b:", min=1, max=Inf, value=1, step=0.1),
				numericInput("kprime_v1", "v1:", min=0, max=Inf, value=1, step=0.001),
				numericInput("kprime_v2", "v2:", min=0, max=Inf, value=2, step=0.001)
			),
								 #2FIX: start here... from lambdap
								 # add warnings about the parameters being recycled ... 
      hr(),
      numericInput("nsamples", "Number of draws:", min = 50, max = 10000, value = 5000, step=50),
      numericInput("randseed", "Rand seed:", min = 1, max = .Machine$integer.max, value = 2015, step=1)
    ,width=3),
    mainPanel(
			tabsetPanel(
				tabPanel("d-d",plotOutput("ddplot")),
				tabPanel("q-q",plotOutput("qqplot")),
				tabPanel("p-p",plotOutput("ppplot"))
			)
		))
,title="Monte Carlo Simulations"))

shinyUI(ui)

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
