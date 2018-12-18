scrapper <- function(URL, CSS){

install.packages('rvest')
library(rvest)

url <- URL
website <- read_html(url)

nodes <- html_nodes(website, CSS)
text <- html_text(nodes)

return(text)

}