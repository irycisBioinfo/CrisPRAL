# serverMosaic_batch

observeEvent(input$go_to_mosaic_batch, {
  
  updateTabsetPanel(session, 'app', selected = 'Mosaic_Finder_batch')
  
})