# app.R  – Stroke Code decision aid
library(shiny)
library(lubridate)

tnk_time_window <- 4.5  # h
tnk_min_nihss   <- 6
evt_min_nihss   <- 6

ui <- fluidPage(
  titlePanel("Stroke Code"),
  sidebarLayout(
    sidebarPanel(width = 4,
      numericInput("age","Age:",value=NA,min=0,max=120),
      radioButtons("sex","Sex:",choices=c("Female","Male"),inline=TRUE),
      textAreaInput("pmh","Past Medical History:"),
      textAreaInput("symptoms","Presenting Symptoms:"),
      textInput("lkw_time","Last Known Well Time (HH:MM AM/PM):",placeholder="02:15 PM"),
      dateInput("lkw_date","Last Known Well Date:"),
      numericInput("nihss","NIHSS Score:",value=NA,min=0,max=42),
      textAreaInput("deficits","Neurologic Deficit(s):"),
      textAreaInput("ct","Head CT Findings:"),
      textAreaInput("cta","CTA Head/Neck Findings:"),
      textAreaInput("ctp","CTP Findings:"),
      checkboxGroupInput("contra","TNK Contraindications:",choices=c(
        "ICH or extensive infarction on NCCT","Ischemic stroke <3 months",
        "Head trauma <3 months","Prior ICH","Intracranial neoplasm",
        "Intracranial/intraspinal surgery <3 months","GI malignancy",
        "GI haemorrhage <3 weeks","Suspected SAH, aortic dissection",
        "SBP >185 / DBP >110 mmHg","Active haemorrhage",
        "Tx-dose LMWH <24 h or DOAC <48 h","Platelets <100 k",
        "INR >1.7 / PT >15 s / PTT >40 s","Glucose <50 mg/dL"))
    ),
    mainPanel(
      h3("Automated Treatment Recommendations"),
      verbatimTextOutput("tnkRec"),
      verbatimTextOutput("evtRec"),
      tags$hr(),
      h4("Rationale"),
      textAreaInput("rationale",NULL,width="100%",height="180px",
                    placeholder="Brief free-text rationale…")
    )
  )
)

server <- function(input, output, session) {
  onset_to_now_min <- reactive({
    req(input$lkw_date,input$lkw_time)
    parsed <- parse_date_time(paste(input$lkw_date,input$lkw_time),
                              orders=c("Y-m-d I:M p","Y-m-d H:M"))
    as.numeric(difftime(Sys.time(), parsed, units="mins"))
  })

  output$tnkRec <- renderText({
    req(input$nihss,input$age)
    if (length(input$contra)>0)
      return("TNK Treatment:  \u2717 Not Recommended  (contraindication present)")
    if (isTRUE(onset_to_now_min()/60 <= tnk_time_window) &&
        input$nihss>=tnk_min_nihss && input$age>=18)
      "TNK Treatment:  \u2713 Recommended"
    else "TNK Treatment:  \u2717 Not Recommended"
  })

  output$evtRec <- renderText({
    req(input$nihss)
    if (input$nihss>=evt_min_nihss)
      "EVT:  \u2713 Recommended  (NIHSS criterion met)"
    else "EVT:  \u2717 Not Recommended"
  })
}

shinyApp(ui, server)
